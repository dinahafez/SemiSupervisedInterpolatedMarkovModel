while read line
do
grep -n “^${line}$” $2;
done < $1
