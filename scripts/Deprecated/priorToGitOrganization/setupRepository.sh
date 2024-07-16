###################
## After inititalizing repo on github
cd ~/Experiments/PhalFNP/
git init
echo '# Compiled source #
###################
*_ignored' > .gitignore
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:NolanBentley/PhalFNP.git
git push -u origin main
