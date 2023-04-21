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
Jeff to show you how if you're interested.

#### Testing

Run `./test.py`

#### Presubmit

Before creating a PR or submitting code, run `./check.sh`.  It will run tests
than check your types and formatting.  This also runs automatically on GitHub
when you make PRs, but it's much faster to catch problems locally first.
