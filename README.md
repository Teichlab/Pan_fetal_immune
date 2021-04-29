# Pan_fetal_immune
Collection of scripts for analysis of pan fetal immune atlas

## Contents

<!-- * `notebooks`: working analysis notebooks
* `utils`: scripts with functions that are used in multiple notebooks
* `manual_annotation`: tables of refined annotations post-integration
* `legacy_code`: old notebooks and scripts that have been refactored in the main analyses, but worth keeping to retain useful snippets

 -->
* [`src`](https://github.com/emdann/Pan_fetal_immune/edit/master/src) contains source code used to analyze the data
* [`metadata`](https://github.com/emdann/Pan_fetal_immune/edit/master/metadata): contains metadata relevant to annotate the samples and pointers to raw and processed data
* `legacy_code`: old notebooks and scripts that have been refactored in the main analyses, but worth keeping to retain useful snippets
* `thrash_n_snippets`: discarded code not included in main analysis

## How to contribute to this repo

To avoid conflicts, each contributor should make a working branch for themselves and then push to master using a pull request:

In your terminal, first clone the repo:
```
git clone emdann/Pan_fetal_immune
```
Then create your personal working branch
```
cd Pan_fetal_immune
git branch working_emma
```
Then switch from master to your working branch
```
git checkout working_emma
```
Now you can add your notebooks and scripts into the relevant folder or add new subfolders, saving your progress committing and pushing to your branch
```
git push origin working_emma
```
When you are happy with adding your progress to the `master` branch you can open a pull request from the github UI, where you can add more detail on what you are adding for all the contributors to read. Then your code will be checked for conflicts and merged to master. 

Of course you can make as many working branches as you like, the important thing is to avoid pushing directly to master as much as possible.

<!-- ### Using issues

Another GitHub feature I find useful is using issues to keep track of what I am working on (example https://github.com/emdann/Pan_fetal_immune/issues/2). When I complete a task I can add a link to the relevant pull request and close the issue. I also use issues to add action items after meetings.

Feel free to do the same, or use issues to discuss ideas and problems (better to keep track than long email threads IMO). 
 -->





