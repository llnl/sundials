---
name: sundials-devs-update-answer-files
description: Update or add SUNDIALS expected-output answer files after intentional test output changes. Use to refresh .out files embedded in examples, update the sundials-codes/answers repository or test/answers submodule, consume GitHub Actions or Jenkins output_files artifacts, configure SUNDIALS_TEST_ANSWER_DIR or SUNDIALS_TEST_OUTPUT_DIR, or follow the documented workflow in doc/developers/testing/Answers.rst.
---

# Update Answer Files

Use this skill when SUNDIALS tests fail because expected-output `.out` files need to change. Treat answer updates as correctness work, not mechanical churn: first confirm the output change is expected and review the numerical/statistical differences before committing any answer file.

Source documentation: `doc/developers/testing/Answers.rst` and `doc/developers/testing/CTest.rst`.

## Workflow

1. Identify the answer-file target.
   - For `.out` files embedded under `examples/`, update the file in the SUNDIALS repository.
   - For platform/compiler answer files, update `sundials-codes/answers`, usually through the `test/answers` submodule.
   - For a new machine directory in `sundials-codes/answers`, use the staged new-machine flow below.

2. Obtain trusted output.
   - Prefer CI-generated output for CI failures because local floating-point differences can be enough to keep CI failing.
   - For GitHub Actions failures, download the `output_files` artifact from the failing workflow or PR check. The unzipped artifact is a SUNDIALS build directory.
   - For Jenkins failures of embedded `examples/*.out`, use the Jenkins build output if local output still fails on Jenkins.
   - For local generation, configure with `-DSUNDIALS_TEST_OUTPUT_DIR=<output-dir>` so test output is written to a predictable directory.

3. Update only the relevant files.
   - Inspect failing test names and paths before copying.
   - Use `scripts/updateOutFiles.py` when updating from a build/output directory because it traverses `examples` and `test/unit_tests` and only updates outputs for failed tests.
   - By default the script reads `<source>/Testing/Temporary/LastTestsFailed.log`; use `--log <dir>` only when the failed-test log is under another build directory.
   - Use `--copy` when intentionally adding a new answer file that does not already exist in the destination. Use `--all` only for an intentionally broad refresh.
   - From the repository root, run:

```bash
cd scripts
./updateOutFiles.py <source-build-or-output-dir> <destination>
```

4. Review diffs before committing.
   - Check every changed `.out` file with `git diff`.
   - Confirm changed solution values, statistics, and tolerances are consistent with the code change.
   - Drop unrelated generated output and avoid broad answer refreshes unless the change intentionally affects many tests.

5. Test with the updated answers.
   - For a local build that should compare against an alternate answer tree, configure with `-DSUNDIALS_TEST_ANSWER_DIR=<answer-dir>`.
   - Run the focused failing tests first, then a broader relevant CTest sweep.

## Embedded Example `.out` Files

Update embedded example answers directly in the SUNDIALS repository.

```bash
cd scripts
./updateOutFiles.py <jenkins-or-local-build-dir> ..
```

If local output passes locally but Jenkins still fails, replace the file with Jenkins-generated output. Getting Jenkins output may require help from the SUNDIALS maintainers named in `doc/developers/testing/Answers.rst`.

## `sundials-codes/answers`

Use this flow for GitHub Actions answer files.

1. Initialize the answer submodule from the SUNDIALS repository root, then enter it.

```bash
git submodule update --init test/answers
cd test/answers
```

If `test/answers` is unavailable or empty, clone `https://github.com/sundials-codes/answers` outside the SUNDIALS source tree instead.

2. Branch from `main`, ideally using the same branch name as the SUNDIALS branch.

```bash
git checkout main
git pull --ff-only
git checkout -b <branch-name>
```

3. Update the relevant precision directory for GitHub Actions:

```text
linux-ubuntu20.04-x86_64/gcc-9.4.0/<double|single|extended>
```

Run the updater from the SUNDIALS repository `scripts` directory:

```bash
cd <sundials-repo>/scripts
./updateOutFiles.py <downloaded-output-files-build-dir> <answers-repo>/linux-ubuntu20.04-x86_64/gcc-9.4.0/<precision>
```

4. If SSH authentication is needed for pushing from the submodule, set:

```bash
git remote set-url origin git@github.com:sundials-codes/answers.git
```

5. Commit the answer updates in the answers repository and open a PR against `sundials-codes/answers`.

## New Machine Answer Directories

Use the staging branch workflow for new platform/compiler directories in `sundials-codes/answers`.

1. Copy the relevant existing answers from `linux-ubuntu20.04-x86_64` into the new machine directory.
2. Commit that initial copy and open a PR to the `staging` branch.
3. After it merges, generate new answers on the new machine, overwrite the copied files, and review the diff.
4. Commit the generated answers and open another PR targeting `staging`.
5. Expect `staging` to be merged into `main` later.

For local generation on `develop`, use a build configured with an explicit output directory:

```bash
git checkout develop
cmake -S . -B build \
  -DSUNDIALS_TEST_ENABLE_DEV_TESTS=ON \
  -DSUNDIALS_TEST_ENABLE_UNIT_TESTS=ON \
  -DSUNDIALS_TEST_OUTPUT_DIR=<machine-output-dir>
cmake --build build -j
ctest --test-dir build --output-on-failure
```

## Gotchas

- If you don't have network access to GitHub, you will need the user to download the `output_files` artifact from the failing workflow or PR check and provide it to you. The unzipped artifact is a SUNDIALS build directory.
