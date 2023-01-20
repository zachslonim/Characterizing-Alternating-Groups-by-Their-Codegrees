# First, we create functions which output the order of a simple group, given its family and parameters
# n is used to represent the rank of group, while q is a prime power representing the order of the base field for lie type groups
function AnOrder(n) # Alternating Groups: n > 4
    return factorial(big(n))/2
end
function AnqOrder(n, q) # Projective Special Linear Groups: n > 0
    n = BigInt(n)
    q = BigInt(q)
    j = 1
    for i in 1:n
        j = j*(q^(i+1)-1)
    end
    return q^((n+1)*n/2)*j/gcd(n+1,q-1)
end
function BnqOrder(n, q) # Orthogonal Groups in Odd Dimension: n > 1
    n = BigInt(n)
    q = BigInt(q)
    j = 1
    for i in 1:n
        j = j*(q^(2*i)-1)
    end
    if q == 2 && n == 2
        return q^(n^2)*j/gcd(2,q-1)/2
    else
        return q^(n^2)*j/gcd(2,q-1)
    end
end
function CnqOrder(n, q) # Symplectic Groups: n > 2
    n = BigInt(n)
    q = BigInt(q)
    j = 1
    for i in 1:n
        j = j*(q^(2*i)-1)
    end
    return q^(n^2)*j/gcd(2,q-1)
end
function DnqOrder(n, q) # Orthogonal Groups in Even Dimension: n > 3
    n = BigInt(n)
    q = BigInt(q)
    j = 1
    for i in 1:n - 1
        j = j*(q^(2*i)-1)
    end
    return (q^n-1)*q^(n*(n-1))*j/gcd(4,q^n-1)
end
function E6qOrder(q) # Chevally Exceptional Group: q is a Prime Power
    q = BigInt(q)
    j = 1
    for i in [2,5,6,8,9,12]
        j = j*(q^i-1)
    end
    return q^36*j/gcd(3,q-1)
end
function E7qOrder(q) # Chevally Exceptional Group: q is a Prime Power
    q = BigInt(q)
    j = 1
    for i in [2,6,8,10,12,14,18]
        j = j*(q^i-1)
    end
    return q^63*j/gcd(2,q-1)
end
function E8qOrder(q) # Chevally Exceptional Group: q is a Prime Power
    q = BigInt(q)
    j = 1
    for i in [2,8,12,14,18,20,24,30]
        j = j*(q^i-1)
    end
    return q^120*j
end
function F4qOrder(q) # Chevally Exceptional Group: q is a Prime Power
    q = BigInt(q)
    j = 1
    for i in [2,6,8,12]
        j = j*(q^i-1)
    end
    return q^24*j
end
function G2qOrder(q) # Chevally Exceptional Group: q is a Prime Power
    q = BigInt(q)
    j = 1
    for i in [2,6]
        j = j*(q^i-1)
    end
    return q^6*j
end
function Anq2Order(n, q) # Projective Special Unitary Groups: n > 1
    n = BigInt(n)
    q = BigInt(q)
    j = 1
    for i in 1:n
        j = j*(q^(i+1)-(-1^(i+1)))
    end
    return q^(n*(n+1)/2)*j/gcd(n+1,q+1)
end
function Dnq2Order(n, q) # Twisted Orthogonal Chevally Groups: n > 3
    n = BigInt(n)
    q = BigInt(q)
    j = 1
    for i in 1:n
        j = j*(q^(2*i)-1)
    end
    return q^(n*(n+1))*(q^n+1)*j/gcd(4,q^n+1)
end
function E6q2Order(q) # Twisted Exceptional Chevally Group: q is a Prime Power
    q = BigInt(q)
    j = 1
    for i in [2,5,6,8,9,12]
        j = j*(q^i-(-1)^i)
    end
    return q^36*j/gcd(3,q+1)
end
function D4q3Order(q)  # Twisted Exceptional Chevally Group: q is a Prime Power
    q = BigInt(q)
    return q^12*(q^8+q^4+1)*(q^6-1)*(q^2-1)
end
function B2q2Order(q) # Suzuki Groups: q is a Prime Power
    q = BigInt(q)
    return q^(2)*(q^(2)+1)*(q-1)
end
function F4q2Order(q) # Ree Groups in Characteristic 2: q = 2 ^ (2m+1)
    q = BigInt(q)
    return q^12*(q^6+1)*(q^4-1)*(q^3+1)*(q-1)
end
function G2q2Order(q) # Ree Groups in Characteristic 3: q = 3 ^ (3m+1)
    q = BigInt(q)
    return q^3*(q^3+1)*(q-1)
end

# Now, we create a function which, given a simple group and an upper bound on its
# number of conjguacy classes, finds which values of n satisfy the order of the simple
# group divides |A_N| and |A_N| is less than the order of the simple group
# multiplied by its number of conjugacy classes
function possibleAn(order, CCNUpperBound)
    possibilities = []
    n = 5
    while true
        if AnOrder(n) > order * CCNUpperBound
            break
        elseif AnOrder(n) % order == 0
            append!(possibilities, n)
        end
        n = n + 1
    end
    return(possibilities)
end

# For the sporadic (or Tits) groups, we can just input each of the 27 orders in turn
# For these, we know the exact number of conjugacy classes from the ATLAS so we use those
println("SPORADIC GROUP POSSIBILITIES")
sporadicOrders = [17971200,7920,95040,443520,10200960,244823040,175560,604800,50232960,86775571046077562880,495766656000,42305421312000,4157776806543360000,64561751654400,4089470473293004800,1255205709190661721292800,44352000,898128000,4030387200,145926144000,448345497600,460815505920,273030912000000,51765179004000000,90745943887872000,4154781481226426191177580544000000,808017424794512875886459904961710757005754368000000000]
sporadicCCNs = [22,10,15,12,17,26,15,21,21,62,42,60,101,65,98,108,24,24,33,36,43,30,54,53,48,184,194]
for i in 1:length(sporadicOrders)
    order = sporadicOrders[i]
    CCN = sporadicCCNs[i]
    possibilities = possibleAn(order, CCN)
    if length(possibilities) > 0
        for n in possibilities
            println(order, " divides |A_", n, "| and |A_", n, "| is less than ", order, "*", CCN, ".")
        end
    else
        println("There are no possible values of n satisfying ", order, " divides |A_n| and |A_n| is less than ", order, "*", CCN, ".")
    end
end
println("")

# For the classical and exceptional groups of Lie type, we first create functions which give upper bounds on the number of conjugacy classes
# These come from tables 1 and 5 in the paper
function AnqCCN(n, q) # Projective Special Linear Groups: n > 0
    n = BigInt(n)
    q = BigInt(q)
    return 2.5 * q^n
end
function BnqCCN(n, q) # Orthogonal Groups in Odd Dimension: n > 1
    n = BigInt(n)
    q = BigInt(q)
    return 7.3 * q^n
end
function CnqCCN(n, q) # Symplectic Groups: n > 2
    n = BigInt(n)
    q = BigInt(q)
    return 15.2 * q^n
end
function DnqCCN(n, q) # Orthogonal Groups in Even Dimension: n > 3
    n = BigInt(n)
    q = BigInt(q)
    return 15 * q^n
end
function E6qCCN(q) # Chevally Exceptional Group: q is a Prime Power
    q = BigInt(q)
    return q^6+q^5+2*q^4+2*q^3+15*q^2+21*q+60
end
function E7qCCN(q) # Chevally Exceptional Group: q is a Prime Power
    q = BigInt(q)
    return q^7+q^7+2*q^5+7*q^4+17*q^3+35*q^2+71*q+103
end
function E8qCCN(q) # Chevally Exceptional Group: q is a Prime Power
    q = BigInt(q)
    return q^8+q^7+2*q^6+3*q^5+10*q^4+16*q^3+40*q^2+67*q+112
end
function F4qCCN(q) # Chevally Exceptional Group: q is a Prime Power
    q = BigInt(q)
    return q^4+2*q^3+7*q^2+15*q+31
end
function G2qCCN(q) # Chevally Exceptional Group: q is a Prime Power
    q = BigInt(q)
    return q^2+2*q+9
end
function Anq2CCN(n, q) # Projective Special Unitary Groups: n > 1
    n = BigInt(n)
    q = BigInt(q)
    return 8.26 * q^n
end
function Dnq2CCN(n, q) # Twisted Orthogonal Chevally Groups: n > 3
    n = BigInt(n)
    q = BigInt(q)
    return 15 * q^n
end
function E6q2CCN(q) # Twisted Exceptional Chevally Group: q is a Prime Power
    q = BigInt(q)
    return q^6+q^5+2*q^4+4*q^3+18*q^2+26*q+62
end
function D4q3CCN(q)  # Twisted Exceptional Chevally Group: q is a Prime Power
    q = BigInt(q)
    return q^4+q^3+q^2+q+6
end
function B2q2CCN(q) # Suzuki Groups: q is a Prime Power
    q = BigInt(q)
    return q+3
end
function F4q2CCN(q) # Ree Groups in Characteristic 2: q = 2 ^ (2m+1)
    q = BigInt(q)
    return q^2+4*q+17
end
function G2q2CCN(q) # Ree Groups in Characteristic 3: q = 3 ^ (3m+1)
    q = BigInt(q)
    return q^2+2*q+9
end

# For the classical groups, once we have the maximum values for m, p, and k, (as explained in the paper)
# we simply loop through all possible combinations and output the ones that are possible exceptions
# We must do this for each of the 6 families of classical lie type groups
# Note that in some cases, the groups output here are not simple but we easily check these and discard before putting them as exceptions

println("CLASSICAL LIE TYPE GROUP POSSIBILITIES")
for m in 1:6 # Projective Special Linear Groups
    for p in [2,3,5,7,11,13,17]
        for k in 1:63
            order = AnqOrder(m, p^k)
            CCN = AnqCCN(m, p^k)
            possibilities = possibleAn(order, CCN)
            if length(possibilities) > 0
                for n in possibilities
                    println("|PSL(", m+1, ",", p^k, ")| divides |A_", n, "| and |A_", n, "| is less than |PSL(", m+1, ",", p^k, ")|*|Irr(PSL(", m+1, ",", p^k, "))|.")
                end
            end
        end
    end
end
for m in 2:2 # Orthogonal Groups in Odd Dimension
    for p in [2,3]
        for k in 1:1
            order = BnqOrder(m, p^k)
            CCN = BnqCCN(m, p^k)
            possibilities = possibleAn(order, CCN)
            if length(possibilities) > 0
                for n in possibilities
                    println("|O(", 2*m+1, ",", p^k, ")| divides |A_", n, "| and |A_", n, "| is less than |O(", 2*m+1, ",", p^k, ")|*|Irr(O(", 2*m+1, ",", p^k, "))|.")
                end
            end
        end
    end
end
for m in 3:4 # Symplectic Groups
    for p in [2]
        for k in 1:2
            order = CnqOrder(m, p^k)
            CCN = CnqCCN(m, p^k)
            possibilities = possibleAn(order, CCN)
            if length(possibilities) > 0
                for n in possibilities
                    println("|PSp(", 2*m, ",", p^k, ")| divides |A_", n, "| and |A_", n, "| is less than |PSp(", 2*m, ",", p^k, ")|*|Irr(PSp(", 2*m, ",", p^k, "))|.")
                end
            end
        end
    end
end
for m in 4:4 # Orthogonal Groups in Even Dimension
    for p in [2]
        for k in 1:1
            order = DnqOrder(m, p^k)
            CCN = DnqCCN(m, p^k)
            possibilities = possibleAn(order, CCN)
            if length(possibilities) > 0
                for n in possibilities
                    println("|O^+(", 2*m, ",", p^k, ")| divides |A_", n, "| and |A_", n, "| is less than |O^+(", 2*m, ",", p^k, ")|*|Irr(O^+(", 2*m, ",", p^k, "))|.")
                end
            end
        end
    end
end
for m in 2:6 # Projective Special Unitary Groups
    for p in [2,3,5,7]
        for k in 1:42
            order = Anq2Order(m, p^k)
            CCN = Anq2CCN(m, p^k)
            possibilities = possibleAn(order, CCN)
            if length(possibilities) > 0
                for n in possibilities
                    println("|PSU(", m+1, ",", p^k, ")| divides |A_", n, "| and |A_", n, "| is less than |PSU(", m+1, ",", p^k, ")|*|Irr(PSU(", m+1, ",", p^k, "))|.")
                end
            end
        end
    end
end
for m in 4:5 # Twisted Orthogonal Chevally Groups
    for p in [2,3]
        for k in 1:3
            order = Dnq2Order(m, p^k)
            CCN = Dnq2CCN(m, p^k)
            possibilities = possibleAn(order, CCN)
            if length(possibilities) > 0
                for n in possibilities
                    println("|O^-(", 2*m, ",", p^k, ")| divides |A_", n, "| and |A_", n, "| is less than |O^-(", 2*m, ",", p^k, ")|*|Irr(O^-(", 2*m, ",", p^k, "))|.")
                end
            end
        end
    end
end
println("")

# For the exceptional lie type groups, we often find that when computing the maximal p and k values,
# we actually get no possible exceptions immediately (i.e. even with q=2^1, the inequality cannot be satisfied).
# As shown in the paper, we have only the following six exceptions to consider here:
# G_2(2)' (which has order half that of G_2(2)), ^3D_4(2), and ^2B_2(2^(2*m+1)) where m = 1, 2, 3, or 4
# For each of these, we just print the number of possible n-values and find it to be 0 in each case:
println("EXCEPTIONAL LIE TYPE GROUP POSSIBILITIES")
println(length(possibleAn(G2qOrder(2)/2, G2qCCN(2))))
println(length(possibleAn(D4q3Order(2), D4q3CCN(2))))
println(length(possibleAn(B2q2Order(2^3), B2q2CCN(2^3))))
println(length(possibleAn(B2q2Order(2^5), B2q2CCN(2^5))))
println(length(possibleAn(B2q2Order(2^7), B2q2CCN(2^7))))
println(length(possibleAn(B2q2Order(2^9), B2q2CCN(2^9))))
