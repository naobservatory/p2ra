## Prevalence to Relative Abundance Project

In this project we're attempting to understand how the prevalence of
human pathogens translates into the relative abundance we see in wastewater
metagenomics.  This work is split across two repositories: in this repo we're
collecting prevalence estimates and building the model, while determining
relative abundances from existing data is in the
[mgs-pipeline](https://github.com/naobservatory/mgs-pipeline) repo.

### Working with prevalence data

In python, run `import pathogens` and then iterate over `pathogens.pathogens`.
Each pathogen implements an `estimate_prevalences` method which gives one or
more estimates.

Run `./summarize.py` to get an overview of the data.

### Statistical model

For an overview of the statistical model see [model.md](model.md).

### Development

#### Making changes

Create a branch named yourname-purpose and push your changes to it.  Then
create a pull request.  Use the "request review" feature to ask for a review:
all PRs need to be reviewed by someone else, and for now include Jeff and Simon
on all PRs unless they're OOO.

Once your change has been approved by your reviewers and passes
[presubmit](#presubmit), you can merge it.  Don't merge someone else's PR
without confirming with them: they may have other changes they've realized they
needed to make, or a tricky branch structure that needs to be resolved in a
particular order.

Handle incoming reviews at least twice a day (morning and afternoon) -- slow
reviews add a lot of friction.  As a PR author you can avoid this friction by
creating another branch that diverges from the code you have under review; ask
Jeff to show you how if you're interested.  Configure [notification
routing](https://github.com/settings/notifications/custom_routing) on github so
that work-related notifications go to your work account.

#### Testing

Run `./test.py`

#### Presubmit

Before creating a PR or submitting code, run `./check.sh`.  It will run tests
than check your types and formatting.  This also runs automatically on GitHub
when you make PRs, but it's much faster to catch problems locally first.

If `./check.sh` complains about formatting or import sorting, you can fix this
automatically with `./check.sh --fix`.

#### Installing pystan

Pystan should be installed along with the other requirements when you run:
```
python -m pip install -r requirements-dev.txt"
```
However, on some non-Linux systems (including M2 Macbooks), one of `pystan`'s dependencies,`httpstan`, may fail to install.
To get around this problem, you can [install httpstan from source](https://httpstan.readthedocs.io/en/latest/installation.html#installation-from-source).
Once it is built and installed, you can then install the requirements file as above.
(Note that you can clone the `httpstan` repo anywhere on your computer.
I recommend doing it outside of the `p2ra` repo directory to that git doesn't try to track it.)
