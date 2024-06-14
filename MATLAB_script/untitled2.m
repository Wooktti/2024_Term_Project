taus = [0,1,4,5,9,10];
D = getDmatrix(taus)

syms t

l0 = (t-1)*(t-4)*(t-5)*(t-9)*(t-10) / ( (0-1)*(0-4)*(0-5)*(0-9)*(0-10) );
l1 = (t-0)*(t-4)*(t-5)*(t-9)*(t-10) / ( (1-0)*(1-4)*(1-5)*(1-9)*(1-10) );
l2 = (t-0)*(t-1)*(t-5)*(t-9)*(t-10) / ( (4-0)*(4-1)*(4-5)*(4-9)*(4-10) );
l3 = (t-0)*(t-1)*(t-4)*(t-9)*(t-10) / ( (5-0)*(5-1)*(5-4)*(5-9)*(5-10) );
l4 = (t-0)*(t-1)*(t-4)*(t-5)*(t-10) / ( (9-0)*(9-1)*(9-4)*(9-5)*(9-10) );
l5 = (t-0)*(t-1)*(t-4)*(t-5)*(t-9) / ( (10-0)*(10-1)*(10-4)*(10-5)*(10-9) );

l0_dot = diff(l0, t);
l1_dot = diff(l1, t);
l2_dot = diff(l2, t);
l3_dot = diff(l3, t);
l4_dot = diff(l4, t);
l5_dot = diff(l5, t);

D = [subs(l0_dot, 1), subs(l1_dot, 1), subs(l2_dot, 1), subs(l3_dot, 1), subs(l4_dot, 1), subs(l5_dot, 1);
    subs(l0_dot, 4), subs(l1_dot, 4), subs(l2_dot, 4), subs(l3_dot, 4), subs(l4_dot, 4), subs(l5_dot, 4);
    subs(l0_dot, 5), subs(l1_dot, 5), subs(l2_dot, 5), subs(l3_dot, 5), subs(l4_dot, 5), subs(l5_dot, 5);
    subs(l0_dot, 9), subs(l1_dot, 9), subs(l2_dot, 9), subs(l3_dot, 9), subs(l4_dot, 9), subs(l5_dot, 9);
    subs(l0_dot, 10), subs(l1_dot, 10), subs(l2_dot, 10), subs(l3_dot, 10), subs(l4_dot, 10), subs(l5_dot, 10)];
D = vpa(D)