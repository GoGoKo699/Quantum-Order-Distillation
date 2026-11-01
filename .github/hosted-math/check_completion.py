#!/usr/bin/env python3
"""Check completion-evidence consistency in an existing private game checkpoint.

No network requests or mutations are made. Evidence is supplied by the launcher;
this checker cannot authenticate observations or certify scientific truth.
See --help and the companion protocol for the evidence contract.
"""
from __future__ import annotations

import argparse
import hashlib
from datetime import datetime, timedelta, timezone
import json
from pathlib import Path, PurePosixPath
from urllib.parse import quote
import re
import sys

CUMULATIVE_STAGES = ("integrated", "release", "cleanup")
STAGES = (*CUMULATIVE_STAGES, "preservation")
SHA256 = re.compile(r"[0-9a-f]{64}\Z")
MAX_REPORT_BYTES = 16 * 1024 * 1024
SHA = re.compile(r"[0-9a-f]{40}\Z")
SAFE = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*\Z")
GAME_ID = re.compile(r"[A-Za-z0-9][A-Za-z0-9_-]{0,127}\Z")
SEMANTIC_TOPICS = {
    "decisions", "results", "failed_approaches", "limitations", "next_steps", "artifacts"
}


def timestamp(value):
    if not isinstance(value, str):
        raise ValueError("timestamp must be a string")
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise ValueError("timestamp needs a timezone")
    return parsed.astimezone(timezone.utc)


class Evidence:
    def __init__(self, now, max_age_hours):
        self.now = now
        self.max_age = timedelta(hours=max_age_hours)

    @staticmethod
    def obj(value, path, errors):
        if not isinstance(value, dict):
            errors.append(f"{path}: required object missing or malformed")
            return {}
        return value

    @staticmethod
    def text(value, path, errors):
        if not isinstance(value, str) or not value.strip():
            errors.append(f"{path}: required nonempty string")
            return ""
        return value

    @staticmethod
    def sha(value, path, errors):
        if not isinstance(value, str) or not SHA.fullmatch(value):
            errors.append(f"{path}: required full lowercase 40-character SHA")
            return ""
        return value

    def observation(self, value, path, errors):
        obj = self.obj(value, path, errors)
        self.text(obj.get("record"), path + ".record", errors)
        try:
            at = timestamp(obj.get("observed_at"))
            if at > self.now + timedelta(minutes=5) or self.now - at > self.max_age:
                errors.append(f"{path}.observed_at: stale or future observation; inspect again")
        except (ValueError, TypeError, OverflowError):
            errors.append(f"{path}.observed_at: required timezone-aware ISO timestamp")
        return obj

    def revision(self, value, path, revision, errors):
        obj = self.observation(value, path, errors)
        for field in ("commit", "tree"):
            actual = self.sha(obj.get(field), path + "." + field, errors)
            matches = (actual in revision.get("allowed_commits", [revision.get("commit")])
                       if field == "commit" else actual == revision.get("tree"))
            if not matches:
                errors.append(f"{path}.{field}: does not match inspected revision")
        return obj

    def merge(self, value, path, errors):
        pr = self.observation(value, path, errors)
        if type(pr.get("number")) is not int or pr["number"] <= 0:
            errors.append(f"{path}.number: required positive PR number")
        if pr.get("merged") is not True:
            errors.append(f"{path}.merged: must be true; a merge hash or closed state is insufficient")
        for field in ("head_sha", "head_tree", "merge_commit_sha", "merge_tree"):
            self.sha(pr.get(field), path + "." + field, errors)
        return pr

    def descendants(self, value, path, tree, errors):
        if not isinstance(value, list):
            errors.append(f"{path}: must be a list (empty is allowed)")
            return []
        edges = []
        for i, raw in enumerate(value):
            key = f"{path}[{i}]"
            edge = self.observation(raw, key, errors)
            source = self.sha(edge.get("from_commit"), key + ".from_commit", errors)
            target = self.sha(edge.get("to_commit"), key + ".to_commit", errors)
            if edge.get("tree") != tree or edge.get("ancestry_confirmed") is not True or edge.get("no_content_change") is not True:
                errors.append(f"{key}: requires observed ancestry and unchanged checked tree")
            else:
                edges.append((source, target))
        return edges

    @staticmethod
    def path_nodes(source, target, edges):
        if (not isinstance(source, str) or not SHA.fullmatch(source)
                or not isinstance(target, str) or not SHA.fullmatch(target)):
            return set()
        if source == target:
            return {source}

        def reachable(start, directed):
            seen = {start}
            while True:
                added = {b for a, b in directed if a in seen} - seen
                if not added:
                    return seen
                seen.update(added)

        forward = reachable(source, edges)
        if target not in forward:
            return set()
        backward = reachable(target, [(b, a) for a, b in edges])
        return forward & backward


def handoff_errors(ev, value, repo, game_id, source, warnings):
    """Validate the same stored handoff for ordinary completion or archival exit."""
    errors = []
    handoff = ev.obj(value, "handoff", errors)
    expected_path = f"projects/{repo.split('/')[-1]}/handoffs/{game_id}/"
    if handoff.get("repository") != "GoGoKo699/Camelot" or handoff.get("private") is not True:
        errors.append("handoff: must identify private GoGoKo699/Camelot")
    if handoff.get("path") != expected_path:
        errors.append("handoff.path: must match this repository and stable game ID")
    path_tree = ev.sha(handoff.get("checked_path_tree"), "handoff.checked_path_tree", errors)
    semantic = ev.observation(handoff.get("semantic"), "handoff.semantic", errors)
    ev.sha(semantic.get("path_tree"), "handoff.semantic.path_tree", errors)
    topics = semantic.get("coverage")
    if not isinstance(topics, list) or any(not isinstance(x, str) for x in topics) or not SEMANTIC_TOPICS.issubset(set(topics)):
        errors.append("handoff.semantic.coverage: missing required important-information topics")
    if semantic.get("complete") is not True or semantic.get("gaps") != []:
        errors.append("handoff.semantic: incomplete or missing gap assessment")
    if (semantic.get("source_commit") != source.get("commit")
            or semantic.get("source_tree") != source.get("tree")
            or semantic.get("path") != expected_path
            or semantic.get("path_tree") != path_tree):
        errors.append("handoff.semantic: assessment must cover this game, checked path tree and final scientific source pin")
    integrity = ev.observation(handoff.get("integrity"), "handoff.integrity", errors)
    if integrity.get("status") != "pass" or integrity.get("scope") != "own_game":
        errors.append("handoff.integrity: passing own-game integrity check required")
    hpr = ev.merge(handoff.get("pull_request"), "handoff.pull_request", errors)
    remote = ev.observation(handoff.get("remote_readback"), "handoff.remote_readback", errors)
    for field in ("commit", "tree", "source_commit", "source_tree"):
        ev.sha(remote.get(field), "handoff.remote_readback." + field, errors)
    if remote.get("path") != expected_path or remote.get("game_id") != game_id or remote.get("content_verified") is not True:
        errors.append("handoff.remote_readback: own-game path, identity and content verification required")
    if (integrity.get("path") != expected_path or integrity.get("path_tree") != path_tree
            or remote.get("path_tree") != path_tree
            or integrity.get("commit") not in (hpr.get("head_sha"), hpr.get("merge_commit_sha"), remote.get("commit"))):
        errors.append("handoff.integrity: checked own-game path tree and revision must match merged remote readback")
    if remote.get("source_commit") != source.get("commit") or remote.get("source_tree") != source.get("tree"):
        errors.append("handoff.remote_readback: stored source pin differs from final scientific main")
    # Shared-store changes are allowed after this PR: verify ancestry and unchanged
    # own-game path instead of requiring the entire Camelot tree to stay frozen.
    if remote.get("commit") != hpr.get("merge_commit_sha"):
        continuation = ev.observation(handoff.get("later_store_revision"), "handoff.later_store_revision", errors)
        if (continuation.get("from_commit") != hpr.get("merge_commit_sha")
                or continuation.get("to_commit") != remote.get("commit")
                or continuation.get("path") != expected_path
                or continuation.get("ancestry_confirmed") is not True
                or continuation.get("own_game_unchanged") is not True):
            errors.append("handoff.later_store_revision: observed ancestry and unchanged own-game path required")
    elif remote.get("tree") != hpr.get("merge_tree"):
        errors.append("handoff.remote_readback.tree: differs from merged store tree")
    if "full_store_audit" in handoff:
        warnings.append("handoff.full_store_audit is informational only; own-game integrity and semantic coverage govern preservation")
    return errors


def safe_relative_path(value):
    return (isinstance(value, str) and bool(value) and "\\" not in value
            and not any(ord(char) < 32 or ord(char) == 127 for char in value)
            and not PurePosixPath(value).is_absolute()
            and all(part not in ("", ".", "..") for part in value.split("/")))


def parse_json(value):
    return json.loads(value, object_pairs_hook=reject_duplicate_keys,
                      parse_constant=lambda x: (_ for _ in ()).throw(ValueError(f"invalid JSON constant: {x}")))


def check_hosted_report(ev, value, repo, revision, required, reports, errors):
    """Bind readiness to inspected report bytes, not a manually copied pass flag."""
    desc = ev.revision(value, "hosted_report", revision, errors)
    path = desc.get("path")
    if not safe_relative_path(path):
        errors.append("hosted_report.path: requires a safe checkpoint-relative file path")
        return
    sha = desc.get("sha256")
    if not isinstance(sha, str) or not SHA256.fullmatch(sha):
        errors.append("hosted_report.sha256: requires complete lowercase SHA256")
    content = reports.get(path) if isinstance(reports, dict) else None
    if not isinstance(content, bytes) or len(content) > MAX_REPORT_BYTES:
        errors.append("hosted_report: actual report bytes are unavailable or exceed 16 MiB")
        return
    if hashlib.sha256(content).hexdigest() != sha:
        errors.append("hosted_report: report bytes do not match sha256")
        return
    try:
        report = parse_json(content)
    except (ValueError, UnicodeError, TypeError) as exc:
        errors.append(f"hosted_report: invalid report JSON: {exc}")
        return
    if not isinstance(report, dict):
        errors.append("hosted_report: report must be a JSON object")
        return
    for field, expected in (("schema_version", 1), ("kind", "hosted_math"),
                            ("repository", repo), ("revision", desc.get("commit")),
                            ("tree", desc.get("tree")), ("renderer_profile", "github"),
                            ("capture_method", "playwright-github"), ("status", "pass"),
                            ("failures", []), ("visual_review_required", True)):
        if (report.get(field) != expected
                or (field == "schema_version" and type(report.get(field)) is not int)
                or (field == "visual_review_required" and report.get(field) is not True)):
            errors.append(f"hosted_report.{field}: invalid or inconsistent collector report")
    ev.observation({"record": path, "observed_at": report.get("captured_at")}, "hosted_report.capture", errors)
    paths = report.get("required_pages")
    if (not isinstance(paths, list) or not paths or not all(safe_relative_path(p) for p in paths)
            or len(paths) != len(set(paths)) or set(paths) != set(required)):
        errors.append("hosted_report.required_pages: must exactly cover declared visual-review pages")
        paths = []
    boundaries = report.get("boundaries")
    if not isinstance(boundaries, list) or any(not isinstance(b, dict) for b in boundaries):
        errors.append("hosted_report.boundaries: needs complete page/source coverage")
        boundaries = []
    seen, count = [], 0
    for boundary in boundaries:
        page = boundary.get("path")
        seen.append(page)
        source_hash = boundary.get("source_sha256")
        if not isinstance(source_hash, str) or not SHA256.fullmatch(source_hash):
            errors.append("hosted_report.boundaries: missing complete source SHA256")
        n, actual = boundary.get("source_expressions"), boundary.get("captured_hosts")
        if type(n) is not int or type(actual) is not int or n < 0 or actual != n:
            errors.append("hosted_report.boundaries: missing or unequal expression/host counts")
        else:
            count += n
        expected_url = f"https://github.com/{repo}/blob/{desc.get('commit')}/{quote(page, safe='/')}" if isinstance(page, str) else None
        if boundary.get("revision") != desc.get("commit") or boundary.get("url") != expected_url:
            errors.append("hosted_report.boundaries: page must identify the exact hosted revision")
    if seen != paths:
        errors.append("hosted_report.boundaries: missing, extra, reordered or duplicate page")
    structure = report.get("structural_comparison")
    if (not isinstance(structure, dict) or structure.get("status") != "pass"
            or type(structure.get("expressions_checked")) is not int
            or structure.get("expressions_checked") != count):
        errors.append("hosted_report.structural_comparison: all source expressions must match")
    views = report.get("viewports")
    if not isinstance(views, list) or len(views) < 2 or any(not isinstance(v, dict) for v in views):
        errors.append("hosted_report.viewports: requires wide and narrow hosted results")
        views = []
    widths = []
    for view in views:
        width, height = view.get("width"), view.get("height")
        if type(width) is not int or type(height) is not int or width <= 0 or height <= 0:
            errors.append("hosted_report.viewports: invalid dimensions")
        else:
            widths.append(width)
        if (view.get("status") != "pass" or type(view.get("structural_matches")) is not int
                or view.get("structural_matches") != count
                or not safe_relative_path(view.get("result_path"))
                or not safe_relative_path(view.get("capture_path"))):
            errors.append("hosted_report.viewports: complete passing reports required at every viewport")
    if not any(w <= 800 for w in widths) or not any(w >= 1200 for w in widths):
        errors.append("hosted_report.viewports: both narrow and wide inspection required")


def preservation_errors(ev, data, repo, game_id, identity_errors, stages, warnings):
    """A user-directed archival exit preserves incomplete work without certifying it."""
    errors = list(identity_errors)
    if data.get("preservation") is None:
        return errors + ["preservation.request: no explicit archival-exit request and evidence supplied"]
    preserve = ev.obj(data.get("preservation"), "preservation", errors)
    request = ev.observation(preserve.get("request"), "preservation.request", errors)
    if (request.get("actor") != "user" or request.get("intent") not in ("preserve_and_stop", "retire_originating_chat")
            or request.get("scope") != "originating_chat" or request.get("game_id") != game_id
            or request.get("repository") != repo):
        errors.append("preservation.request: actual user request for this game's originating chat required")
    source = ev.observation(preserve.get("source"), "preservation.source", errors)
    for field in ("commit", "tree"):
        ev.sha(source.get(field), "preservation.source." + field, errors)
    if source.get("repository") != repo:
        errors.append("preservation.source: repository must match this game")
    integration = data.get("integration")
    if isinstance(integration, dict) and "main" in integration:
        main = ev.obj(integration["main"], "integration.main", errors)
        if any(main.get(k) != source.get(k) for k in ("commit", "tree")):
            errors.append("preservation.source: conflicts with observed scientific main")
    errors.extend(handoff_errors(ev, data.get("handoff"), repo, game_id, source, warnings))
    handoff = data.get("handoff") if isinstance(data.get("handoff"), dict) else {}
    resume = ev.observation(preserve.get("resumption"), "preservation.resumption", errors)
    for field, expected in (("game_id", game_id), ("repository", repo),
                            ("source_commit", source.get("commit")), ("source_tree", source.get("tree")),
                            ("path", handoff.get("path")), ("path_tree", handoff.get("checked_path_tree"))):
        if resume.get(field) != expected:
            errors.append(f"preservation.resumption.{field}: must match the saved game and source")
    for flag, expected in (("self_contained", True), ("requires_original_chat", False),
                           ("artifacts_reconciled", True), ("unfinished_work_reviewed", True)):
        if resume.get(flag) is not expected:
            errors.append(f"preservation.resumption.{flag}: required {expected}")
    ev.text(resume.get("artifact_record"), "preservation.resumption.artifact_record", errors)
    if resume.get("missing_essential_artifacts") != []:
        errors.append("preservation.resumption: essential artifacts remain missing or unassessed")
    tasks = resume.get("unfinished_tasks")
    if not isinstance(tasks, list) or any(not isinstance(t, dict) for t in tasks):
        errors.append("preservation.resumption.unfinished_tasks: requires an explicit list")
        tasks = []
    task_stages = set()
    for i, task in enumerate(tasks):
        label = f"preservation.resumption.unfinished_tasks[{i}]"
        if task.get("stage") not in CUMULATIVE_STAGES:
            errors.append(label + ".stage: unknown remaining stage")
        else:
            task_stages.add(task["stage"])
        for key in ("reason", "owner", "next_action"):
            ev.text(task.get(key), label + "." + key, errors)
    # A failed scientific or display review can be preserved, not silently erased.
    if not stages["integrated"]["ok"] and "integrated" not in task_stages:
        errors.append("preservation.resumption: unfinished integration/science needs a concrete continuation task")
    if not stages["release"]["ok"] and "release" not in task_stages:
        errors.append("preservation.resumption: unfinished release review needs a concrete continuation task")
    block = data.get("blocker")
    if isinstance(block, dict):
        if block.get("stage") == "cleanup":
            errors.append("preservation: own-game storage blocker remains unresolved")
        elif block.get("stage") in ("integrated", "release"):
            if not any(all(task.get(k) == block.get(k) for k in ("stage", "reason", "owner", "next_action")) for task in tasks):
                errors.append("preservation.resumption: retain the exact unresolved blocker and permitted next action")
        else:
            errors.append("preservation: malformed blocker cannot be discarded")
    elif block is not None:
        errors.append("preservation: malformed blocker cannot be discarded")
    return errors


def evaluate(checkpoint, *, now=None, max_age_hours=24, reports=None):
    """Return cumulative completion and independent preservation evidence verdicts."""
    if now is None:
        now = datetime.now(timezone.utc)
    if now.tzinfo is None or not isinstance(max_age_hours, (int, float)) or isinstance(max_age_hours, bool) or not 0 < max_age_hours <= 168:
        raise ValueError("now needs a timezone; max_age_hours must be in (0, 168]")
    ev = Evidence(now.astimezone(timezone.utc), max_age_hours)
    base, release, cleanup, warnings = [], [], [], []
    root = ev.obj(checkpoint, "checkpoint", base)
    data = ev.obj(root.get("completion_evidence"), "completion_evidence", base)
    if type(data.get("schema_version")) is not int or data["schema_version"] != 1:
        base.append("completion_evidence.schema_version: must be 1")
    game = ev.obj(data.get("game"), "game", base)
    game_id = ev.text(game.get("id"), "game.id", base)
    repo = ev.text(game.get("repository"), "game.repository", base)
    if not GAME_ID.fullmatch(game_id):
        base.append("game.id: must match Camelot's 1-128 character alphanumeric/underscore/hyphen identifier")
    if len(repo.split("/")) != 2 or not all(SAFE.fullmatch(x) for x in repo.split("/")):
        base.append("game.repository: must be owner/name")
    ev.text(game.get("grant_record"), "game.grant_record", base)
    for outer, inner in (("game_id", game_id), ("repository", repo)):
        if outer in root and root[outer] != inner:
            base.append(f"game: {outer} disagrees with existing checkpoint")
    identity_errors = list(base)
    revision = ev.obj(data.get("revision"), "revision", base)
    for field in ("commit", "tree"):
        ev.sha(revision.get(field), "revision." + field, base)
    integration = ev.obj(data.get("integration"), "integration", base)
    pr = ev.merge(integration.get("pull_request"), "integration.pull_request", base)
    main = ev.observation(integration.get("main"), "integration.main", base)
    for field in ("commit", "tree"):
        ev.sha(main.get(field), "integration.main." + field, base)
    edges = ev.descendants(integration.get("unchanged_descendants", []), "integration.unchanged_descendants", revision.get("tree"), base)
    allowed_commits = set()
    for label, source, target in (
        ("checked revision to PR head", revision.get("commit"), pr.get("head_sha")),
        ("merged PR to current main", pr.get("merge_commit_sha"), main.get("commit")),
    ):
        nodes = ev.path_nodes(source, target, edges)
        if not nodes:
            base.append(f"integration: {label} lacks an observed no-content-change descendant path")
        allowed_commits.update(nodes)
    if any(x != revision.get("tree") for x in (pr.get("head_tree"), pr.get("merge_tree"), main.get("tree"))):
        base.append("integration: checked, PR-head, merged and current-main trees must match")
    # The PR head, actual merge and documented unchanged descendants can all be
    # inspected separately. Their shared tree is what binds checks and captures.
    checked_revision = {**revision, "allowed_commits": list(allowed_commits)}
    policy = ev.observation(data.get("repository_policy"), "repository_policy", base)
    if type(policy.get("ci_required")) is not bool:
        base.append("repository_policy.ci_required: required boolean from inspected policy")
    checks = data.get("checks")
    if not isinstance(checks, list):
        base.append("checks: required list")
        checks = []
    kinds = {"science": base, "ci": base, "learning": release, "hosted": release}
    found = set()
    for i, raw in enumerate(checks):
        path = f"checks[{i}]"
        if not isinstance(raw, dict) or not isinstance(raw.get("kind"), str) or raw["kind"] not in kinds:
            base.append(f"{path}: unknown or missing check kind")
            continue
        kind = raw["kind"]
        found.add(kind)
        dest = kinds[kind]
        check = ev.revision(raw, path, checked_revision, dest)
        ev.text(check.get("name"), path + ".name", dest)
        status = check.get("status")
        if status not in ("pass", "fail", "blocked"):
            dest.append(f"{path}.status: must be pass, fail or blocked")
        elif status != "pass":
            ev.text(check.get("reason"), path + ".reason", dest)
            if kind == "ci" and status == "blocked" and policy.get("ci_required") is False:
                warnings.append(f"{path}: optional CI remains {status}; local checks do not change that status")
            else:
                dest.append(f"{path}: required {kind} check is {status}")
    for kind in ("science", "learning", "hosted"):
        if kind not in found:
            kinds[kind].append(f"checks: missing required {kind} evidence")
    if policy.get("ci_required") is True and "ci" not in found:
        base.append("checks: repository policy requires CI evidence")
    ev.text(game.get("anchor"), "game.anchor", release)
    visual = ev.obj(data.get("visual_review"), "visual_review", release)
    required = visual.get("required_pages")
    if (not isinstance(required, list) or not required or not all(safe_relative_path(x) for x in required)
            or len(required) != len(set(required))):
        release.append("visual_review.required_pages: required nonempty unique safe page-path list")
        required = []
    ev.text(visual.get("scope_record"), "visual_review.scope_record", release)
    pages = visual.get("pages")
    if not isinstance(pages, list):
        release.append("visual_review.pages: required list of actual captures and source comparisons")
        pages = []
    seen = set()
    for i, raw in enumerate(pages):
        path = f"visual_review.pages[{i}]"
        page = ev.revision(raw, path, checked_revision, release)
        name = ev.text(page.get("path"), path + ".path", release)
        if not safe_relative_path(name):
            release.append(f"{path}.path: requires a safe repository-relative page path")
        if name not in required:
            release.append(f"{path}.path: page is outside declared hosted/visual-review scope")
        if name in seen:
            release.append(f"{path}: duplicate page")
        seen.add(name)
        if page.get("source_compared") is not True or page.get("display_inspected") is not True:
            release.append(f"{path}: source comparison and actual visual inspection both required")
    for path in required:
        if path not in seen:
            release.append(f"visual_review: missing required page {path}")
    check_hosted_report(ev, data.get("hosted_report"), repo, checked_revision, required, reports or {}, release)
    cleanup.extend(handoff_errors(ev, data.get("handoff"), repo, game_id, main, warnings))
    blocker = data.get("blocker")
    if blocker is not None:
        block = ev.obj(blocker, "blocker", base)
        target = block.get("stage")
        if target not in CUMULATIVE_STAGES:
            base.append("blocker.stage: must name integrated, release or cleanup")
        else:
            errors = {"integrated": base, "release": release, "cleanup": cleanup}[target]
            for key in ("reason", "owner", "next_action"):
                ev.text(block.get(key), "blocker." + key, errors)
            errors.append(f"blocker: {block.get('reason', 'unresolved')} (owner: {block.get('owner', 'missing')}; next: {block.get('next_action', 'missing')})")
    stages = {}
    cumulative = []
    for name, errors in zip(CUMULATIVE_STAGES, (base, release, cleanup)):
        cumulative += errors
        stages[name] = {"ok": not cumulative, "errors": list(dict.fromkeys(cumulative))}
    preserve = preservation_errors(ev, data, repo, game_id, identity_errors, stages, warnings)
    stages["preservation"] = {"ok": not preserve, "errors": list(dict.fromkeys(preserve))}
    return {"schema_version": 1, "preservation_ready": not preserve,
            "science_integrated": stages["integrated"]["ok"],
            "release_ready": stages["release"]["ok"], "chat_cleanup_ready": stages["cleanup"]["ok"],
            "stages": stages, "warnings": warnings,
            "scope": "Evidence consistency only; no authenticity, scientific-truth or exhaustive-coverage certification."}


def reject_duplicate_keys(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def load_report_bytes(checkpoint_path, checkpoint):
    """Read only a bounded regular report inside the checkpoint's own directory."""
    data = checkpoint.get("completion_evidence") if isinstance(checkpoint, dict) else None
    desc = data.get("hosted_report") if isinstance(data, dict) else None
    path = desc.get("path") if isinstance(desc, dict) else None
    if not safe_relative_path(path):
        return {}
    try:
        base = checkpoint_path.resolve().parent
        candidate = (base / path).resolve()
        if not candidate.is_relative_to(base):
            return {}
        if not candidate.is_file() or candidate.stat().st_size > MAX_REPORT_BYTES:
            return {}
        with candidate.open("rb") as stream:
            content = stream.read(MAX_REPORT_BYTES + 1)
        return {path: content} if len(content) <= MAX_REPORT_BYTES else {}
    except (OSError, RuntimeError, ValueError):
        return {}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("checkpoint", type=Path, help="Existing private checkpoint containing completion_evidence")
    parser.add_argument("--stage", choices=STAGES, required=True)
    parser.add_argument("--now", help="Timezone-aware ISO timestamp (default: current UTC)")
    parser.add_argument("--max-age-hours", type=float, default=24,
                        help="Maximum age of supplied observations, at most 168 hours (default 24)")
    args = parser.parse_args()
    try:
        data = parse_json(args.checkpoint.read_text(encoding="utf-8"))
        reports = load_report_bytes(args.checkpoint, data)
        result = evaluate(data, now=timestamp(args.now) if args.now else None,
                          max_age_hours=args.max_age_hours, reports=reports)
    except (OSError, UnicodeError, ValueError, TypeError, OverflowError) as exc:
        print(json.dumps({"ok": False, "error": str(exc)}))
        return 2
    print(json.dumps(result, indent=2))
    return 0 if result["stages"][args.stage]["ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
