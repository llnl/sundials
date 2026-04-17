---
name: devs-address-pr-comments
description: Fetch and apply GitHub pull-request review comments locally, especially when the task includes reconciling new review traffic, following reply chains from specific reviewers/authors, and keeping a scratch notes file that maps comment IDs to fixes.
---

# Address PR review comments

Use this skill when the user asks to pull review comments from a GitHub PR, apply the requested changes in the local branch, and track what was fixed.

Prefer this workflow over ad hoc browsing when any of these are true:

- The user gives a PR URL or repo/PR number and wants comments applied locally.
- There has been new PR traffic and the branch already contains earlier comment-driven edits.
- The user wants special attention paid to replies from a specific author.
- The user wants a notes file in `.agents/scratch/` mapping comment IDs to fixes.

## Inputs to pin down

- Repo and PR number.
- Whether to fetch only inline review comments or also general issue comments.
- Any author whose replies should override earlier suggestions.
- Notes file path. Default to `.agents/scratch/pr-<PR>-review-notes.md`.

If the user already supplied these, do not ask again.

## Workflow

1. Refresh the PR review comments into scratch storage.

   Prefer the GitHub REST API because `gh` may not be installed:

   ```bash
   curl -fsSL -H 'Accept: application/vnd.github+json' \
     "https://api.github.com/repos/<owner>/<repo>/pulls/<pr>/comments?per_page=100" \
     > .agents/scratch/pr-<pr>-comments-latest.json
   ```

   If `gh` is available, `gh api repos/<owner>/<repo>/pulls/<pr>/comments --paginate` is fine.
   If `gh` is not available, suggest that the user installs it.

2. Extract the actionable threads.

   Use `jq` to inspect:

   - `id`
   - `in_reply_to_id`
   - `user.login`
   - `path`
   - `body`
   - `diff_hunk`

   Pay attention to reply chains. A later reply from the PR author or an explicitly named person may supersede the original suggestion.

3. Reconcile with the local branch state before editing.

   - Read the files named in the comment threads.
   - Read any existing scratch notes and current `git status`.
   - Do not assume a comment is still outstanding just because it appears in the PR.
   - If an earlier local edit conflicts with a newer reply, follow the newer decision and update the notes.

4. Apply the changes locally.

   - Use `apply_patch` for manual edits.
   - Do not revert unrelated user changes.
   - Favor the narrowest change that resolves the reviewed concern.
   - When comments describe naming/API cleanups, update call sites, docs, and bindings together so the branch stays internally consistent.
   - Make a commit (but do not push) for every comment that you address. The commit message format should be:
   
   ```text
   Agent (<your name and version>) generated commit:

   <URL of comment(s) being addressed>
   
   <what was done and why>
   ```

   Be succinct in your explanation.

5. Maintain the notes file while working.

   Record one flat bullet per addressed thread, including:

   - comment ID or IDs
   - short summary of the requested change
   - brief note on what was changed locally

   If a later reply reversed an earlier change, record both IDs and explain the final decision.

6. Verify the result.

   Minimum checks:

   - `git diff --check`
   - targeted build or tests for touched areas

   In SUNDIALS, prefer targeted developer builds in an out-of-source build dir such as:

   ```bash
   cmake -S . -B build-pr-review \
     -DCMAKE_BUILD_TYPE=Debug \
     -DBUILD_TESTING=ON \
     -DSUNDIALS_TEST_ENABLE_UNIT_TESTS=ON \
     -DSUNDIALS_TEST_ENABLE_DEV_TESTS=ON
   cmake --build build-pr-review -j 4 --target <relevant-targets>
   ```

   Run a focused example or test binary when the touched code has one.

## Review heuristics

- Treat direct replies from the branch author or explicitly prioritized reviewer as tie-breakers.
- Prefer the latest concrete design decision over earlier broad suggestions.
- Watch for comments that change semantics, not just style.

## Output expectations

Before finishing, ensure:

- the local code reflects the latest actionable review guidance,
- `.agents/scratch/pr-<pr>-comments-latest.json` exists or was refreshed,
- the notes file in `.agents/scratch/` is updated,
- if GitHub replies were posted, each one clearly says it was agent-generated,
- the final response summarizes what changed and what verification ran.
