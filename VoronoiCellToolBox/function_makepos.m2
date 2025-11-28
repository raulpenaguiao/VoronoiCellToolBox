-- it takes a polynomial, a list of values, and flips the sign of the polynomial if
-- in the evaluation in these values it is negative
-- makepos(polynom, lvalues, d, RingR)
-- Input:
--   polynom: a polynomial (rational function)
--   lvalues: list of values for substitution
--   d: dimension
--   RingR: polynomial ring
-- Output:
--   The polynomial, possibly with sign flipped to ensure positivity at lvalues
-- Description:
--   Substitutes lvalues into polynom and flips sign if result is negative.
-- Example:
--   makepos(q_0+q_1, {2,-1}, 2, QQ[q_0,q_1])
--   -- Output: q_0+q_1 or -(q_0+q_1) depending on sign
makepos = (polynom, lvalues, d, RingR) -> (
    G = d*(d+1)//2;
    slst = {};
    genList = gens RingR;
    for i from 0 to G-1 do (
        slst = append(slst, genList_i => lvalues_i);
    );
    subvalue1 = sub(numerator polynom, slst);
    subvalue2 = sub(denominator polynom, slst);
    subvalue = subvalue1/subvalue2;
    if subvalue < 0 then return -polynom
    else return polynom
);
