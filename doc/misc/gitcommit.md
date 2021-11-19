git checkout -b bname
git commit -a -m 'msgs'
git push origin bname
git checkout master
git fetch upstream
git merge upstream/master
git push
git branch -D bname
git push origin --delete bname
