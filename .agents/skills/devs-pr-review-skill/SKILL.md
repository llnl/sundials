---
name: devs-pr-review-skill
description: Review GitHub pull requests in a local checkout, draft actionable review comments, and optionally publish selected comments with `gh`. Use when the user asks to review a PR, inspect a GitHub diff against its base branch, prepare markdown review notes before posting, convert findings into GitHub review comments, or publish only a chosen subset of drafted comments while identifying the agent in each posted comment.
---

# Devs PR Review Skill

## Overview

Review the PR against its actual base branch, focus first on bugs, regressions, missing tests, and API/documentation mismatches, and keep the output split between drafted comments and explicitly approved published comments.

Prefer concrete review artifacts over vague summaries: file paths, affected behavior, missing coverage, and the exact comment text that should be posted.

## Workflow

### 1. Establish PR Context

Use `gh` and local git state before forming conclusions.

- Check the current branch and worktree status.
- Inspect remotes so the PR's upstream repo is unambiguous.
- Read PR metadata with `gh pr view <number> --repo <owner>/<repo> --json ...`.
- Fetch the PR head and the base branch locally.
- Compare `base...pr-head`, not the current checkout against itself.

Typical commands:

```bash
gh pr view 916 --repo LLNL/sundials --json number,title,headRefName,baseRefName,files,commits,url
git fetch origin pull/916/head:pr-916
git fetch origin develop
git diff --stat origin/develop...pr-916
```

### 2. Review the Highest-Risk Changes First

Prioritize:

- Core runtime logic and API changes
- Lifetime and ownership changes
- Cross-language bindings and generated interfaces
- Tests and documentation coverage for the new behavior

Look for:

- Behavior that contradicts the new API contract
- Missing reset/unset paths
- Missing negative-path or regression coverage
- Binding changes without corresponding language-level tests
- Scope creep unrelated to the main feature

If practical, run targeted validation. If the environment is missing required tools, state that plainly in the draft.

### 3. Draft Review Comments Before Posting

When the user asks to prepare comments but not submit them yet, create a markdown draft file in the current workspace.

Structure the draft with:

- PR identifier and status
- Short review summary
- Blocking comments
- Suggestion comments
- Submission plan

For each drafted comment include:

- Target file
- Severity
- Draft comment text
- Short rationale for why it should be posted

Keep comments actionable and specific. Prefer one issue per comment.

### 4. Publish Only User-Approved Comments

Do not publish review comments unless the user explicitly asks.

When publishing:

- Post only the numbered or named comments the user approved.
- Omit any comments the user said not to publish.
- Identify the agent in each posted comment.
- Prefer inline review comments when the API supports the requested positioning.
- If the endpoint rejects modern line-based positioning, fall back to file-level PR comments rather than abandoning the comment entirely.

Agent signature format:

```text
Agent: Codex (GPT-5)
```

If the user asks for a different name/version string, follow that instead.

### 5. Use `gh` Pragmatically

Useful patterns:

```bash
gh pr view <number> --repo <owner>/<repo> --json ...
gh pr diff <number> --repo <owner>/<repo> --patch
gh api repos/<owner>/<repo>/pulls/<number>/comments -X POST ...
```

For file-level comments, the PR comments API accepts `subject_type=file`.

If inline review-comment positioning fails, do not loop blindly. Read the error, adjust once, and fall back to a file-level comment if needed.

## Review Standard

Default to a code-review mindset:

- Findings first
- Ordered by severity
- Focus on correctness, regressions, and missing tests
- Keep summaries brief

Good review comments explain:

- What is wrong
- Why it matters
- What change would resolve it

Avoid:

- Generic praise comments
- Restating the diff without analysis
- Posting speculative comments without a concrete failure mode or test gap

## Output Rules

When the user asks for a draft only:

- Create a markdown file
- Do not submit anything to GitHub
- Tell the user where the file is

When the user asks to publish:

- Report which comments were posted
- Link to the resulting PR discussion URLs when available
- State clearly if any draft comments were intentionally not published
