cd "~/Dropbox/Alia/reg_copy"
do "reg_copy.ado"

sysuse auto
reg_copy price mpg trunk displacement if rep78 == 3
reg price mpg trunk displacement if rep78 == 3

