spawn ./tecout.ex

set timeout 1

expect {
    "=====  CONVERT FDNS OUTPUT TO TECPLOT OUTPUT =====" {exp_continue}
    " ==> SELECT THE FILE TYPE OF FDNS GRID/FLOW FILES" {exp_continue}
    "     0: BINARY;  1: FORMATTED" {send "0\r"}
}

expect "==> INPUT FILENAME FOR FDNS INPUT CARD"
send "fort.11\r"

expect "==> INPUT FILENAME FOR FDNS GRID FILE"
send "fort.12\r"

expect "==> INPUT FILENAME FOR FDNS FLOW FILE"
send "fort.23\r"

expect "==> INPUT FILENAME FOR TECPLOT OUTPUT FILE"
send "tecout.dat\r"

expect " ==> SELECT THE INTERVAL FOR TECPLOT OUTPUT(3-D only)"
send "1 1 1\r"

expect eof

