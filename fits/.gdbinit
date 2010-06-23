set print pretty on
set args 'dvcs' 'thi' 'thi'

# macro for jumping
define bj
    tbreak $arg0
    jump $arg0
end

define prkin
    x/3fg kinematics_
end
define prcff
    x/3fg cff_
end
define prpar
    x/30fg par_
end

#break MAIN__
#break GepardInitInternal 
#break  F2Internal
#break  cffHInternal
break cfff_
#break fcn_
#break evolc_
#break BCAInternal
#run -linkmode listen -linkname 5000
run
