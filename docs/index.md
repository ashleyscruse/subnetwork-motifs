---
title: Getting Started
layout: default
---

<div class="page-nav">
  <a href="./" class="current">Getting Started</a>
  <a href="scope.html">Your Scope</a>
</div>

# Getting Started

Welcome to the subnetwork-motifs project. This page walks you through setup end to end. Follow the steps in order. You should reach the "Verify" step in about 15 to 25 minutes on a laptop, most of it waiting for `pip install` and `torch` to download.

## Prerequisites

Before you start, make sure you have:

- A GitHub account (send me your username so I can add you as a collaborator).
- A terminal you are comfortable using.
- `git` installed. Check with `git --version`.
- Python 3.13. Check with `python3.13 --version`. If it is missing:
  - macOS: `brew install python@3.13`
  - Ubuntu / Debian: `sudo apt install python3.13 python3.13-venv`

## Step 1: Accept the invitation

Check your email for a GitHub invitation to the `subnetwork-motifs` repository and click "Accept invitation." You will not be able to clone or push until you accept it.

## Step 2: Set up SSH (recommended)

Long term, SSH is easier than typing a token every time you push. Skip to Step 3 if you already have an SSH key on your GitHub account.

Check whether you already have a key:

```bash
ls -la ~/.ssh
```

If you do not see `id_ed25519` or `id_rsa`, generate one:

```bash
ssh-keygen -t ed25519 -C "your-email@example.com"
```

Press Enter to accept the defaults. Optionally set a passphrase.

Copy the public key:

```bash
cat ~/.ssh/id_ed25519.pub
```

On GitHub, go to Settings, then SSH and GPG keys, then New SSH key, and paste the key. Verify the connection:

```bash
ssh -T git@github.com
```

You should see a message like "Hi username! You've successfully authenticated..."

If you prefer HTTPS, skip this step and generate a personal access token instead at [https://github.com/settings/tokens](https://github.com/settings/tokens) with the `repo` scope. GitHub will prompt for the token the first time you push.

## Step 3: Clone the repo

Choose where you want the project to live:

```bash
cd ~/Documents
git clone git@github.com:ashleyscruse/subnetwork-motifs.git
cd subnetwork-motifs
```

If you are using HTTPS instead of SSH:

```bash
git clone https://github.com/ashleyscruse/subnetwork-motifs.git
```

## Step 4: Switch to your working branch

You have your own branch, `delva`. Do all your work there.

```bash
git checkout delva
```

You should see a message that the branch is now tracking `origin/delva`.

## Step 5: Create a virtual environment and install dependencies

```bash
python3.13 -m venv venv
source venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

The last line takes a few minutes because `torch` is large. When it finishes, your prompt should show `(venv)` at the beginning.

## Step 6: Reproduce the dataset

Run each script in order. The data itself is not committed to the repo, only the scripts that reproduce it.

```bash
python -m src.data.download_harbison
python -m src.data.parse_harbison
python -m src.data.download_wapinski
python -m src.data.parse_wapinski
python -m src.data.download_gasch
python -m src.data.parse_gasch --long
python -m src.data.download_yeast_proteome
python -m src.data.compute_paralog_similarity
python -m src.data.build_graph
```

Each script prints a summary at the end. Total time is about 10 minutes.

## Step 7: Verify

Load the consolidated graph artifact and check its shape:

```bash
python -c "import torch; d = torch.load('data/processed/yeast_graph.pt', weights_only=False); print(f'nodes: {d.num_nodes}, edges: {d.edge_index.shape[1]}, pairs: {d.paralog_edges.shape[1]}')"
```

Expected output:

```
nodes: 6422, edges: 10853, pairs: 801
```

If you see those numbers, you are set up.

## GitHub workflow

The rules for working on this repo:

1. Never push directly to `main`. Always work on your `delva` branch.
2. Commit often on `delva`, with short, clear messages.
3. When you want changes merged, open a Pull Request from `delva` into `main` on GitHub. I review before merging.
4. Do not commit secrets: no passwords, no personal information, no API keys.
5. Read `DELVA.md` and `CONTRIBUTING.md` in the repo for reference.
6. Ask questions early. Do not sit on a blocker.

## Making your first commit

Try a small change end to end so you know the workflow. For example, add a note to one of the notebooks or edit `DELVA.md`.

```bash
git status                            # see what changed
git add <the-file>                    # stage it
git commit -m "First commit: verified environment"
git push origin delva                 # push to your branch
```

After you push, go to the repo on GitHub. You will see a "Compare & pull request" button at the top of the page. Click it, add a short description, and submit the PR. I will review.

## Working with AI assistants

If you use Claude, Cursor, Antigravity, or a similar tool:

- Ask the assistant to read `docs/plan.md` and `docs/data.md` first, so it has real project context before suggesting changes.
- Keep it on the `delva` branch.
- Do not let it push directly to `main` or open pull requests on your behalf without your review.
- Verify code suggestions against the actual scripts before running them.

## Troubleshooting

**Torch install fails on Apple Silicon:**

```bash
pip install torch --index-url https://download.pytorch.org/whl/cpu
```

**Torch-geometric errors during install:**

```bash
pip install torch-geometric --no-cache-dir
```

**Torch complains about Python 3.14:** use `python3.13` explicitly when creating the venv.

**SSH permission denied when cloning:** verify your SSH key is added to GitHub and retest with `ssh -T git@github.com`.

**`venv/bin/activate: no such file`:** you probably ran `python3.13 -m venv venv` from the wrong directory. Confirm you are inside the `subnetwork-motifs/` folder.

## Where to look next

- `docs/plan.md`: full project plan and phase list.
- `docs/data.md`: data documentation, sources, and schemas.
- `CONTRIBUTING.md`: environment setup, workflow, code style.
- `DELVA.md`: your onboarding note at the repo root.
- `results/predicted_pi.csv`: baseline predictions from the first end-to-end training run.

Questions? Message me.
