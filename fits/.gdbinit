set print pretty on
#set args 0

# macro for jumping
define bj
    tbreak $arg0
    jump $arg0
end

define prparint
    x/20d parint_
end
define prcff
    x/6f cff_
end

break MAIN__
break  F2Internal
break  cffHInternal
run -linkmode listen -linkname 5000
