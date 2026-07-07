# Onboarding: Teschon Delva

Welcome to the subnetwork-motifs project.

## Start here

Read [`CONTRIBUTING.md`](CONTRIBUTING.md) for the full setup and workflow. In brief:

1. Clone the repo and check out your branch:
   ```bash
   git clone git@github.com:ashleyscruse/subnetwork-motifs.git
   cd subnetwork-motifs
   git checkout delva
   ```
2. Follow `CONTRIBUTING.md` sections 2 through 4 to set up the environment and reproduce the data.

## Your branch

You have your own branch, `delva`. Work there. Do not push directly to `main`. When you want changes merged, open a pull request from `delva` into `main` on GitHub and I will review it.

## Where to look

- `docs/plan.md`: full project plan and phase list.
- `docs/data.md`: data documentation, sources, and schemas.
- `CONTRIBUTING.md`: environment setup, workflow, code style.
- `results/predicted_pi.csv`: baseline predictions from the first end-to-end training run.

## Working with AI assistants

If you use Claude, Cursor, Antigravity, or similar tools:

- Ask the assistant to read `docs/plan.md` and `docs/data.md` first, so it has real project context before suggesting changes.
- Keep it on the `delva` branch.
- Do not let it push directly to `main` or open pull requests on your behalf without your review.
- Verify code suggestions against the actual scripts before running them.

## Questions

Message me and I will get back to you.
