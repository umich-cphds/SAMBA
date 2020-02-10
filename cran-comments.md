This is a re-re-submission.

## Test environments
* Ubuntu 16.04 (travis): r-devel (2020-01-28 r77738), r-release (3.6.2),
    r-oldrelease (3.5.3)
* winbuilder: r-devel (2020-01-28 r77738), r-release (3.6.2)

## R CMD check results
> There were no ERRORs or WARNINGs.

> There was 1 NOTE:

> Maintainer: ‘Alexander Rix <alexrix@umich.edu>’

This is a new submission.

## Reviewer comments

> Please always write package names, software names and API names in
> single quotes in title and description. e.g: --> 'SAMBA'

Fixed

> You are setting options(warn=-1) in your function/vignette/example. This
> is not allowed.

Fixed

> \dontrun{} should only be used if the example really cannot be executed
> (e.g. because of missing additional software, missing API keys, ...) by
> the user. That's why wrapping examples in \dontrun{} adds the comment
> ("# Not run:") as a warning for the user.
> Does not seem necessary.
> Please replace \dontrun with \donttest.

Fixed

## Previous Reviewer comments

> Found the following (possibly) invalid file URI:
>   URI: commits/master
>     From: README.md

I believe this has been fixed.
