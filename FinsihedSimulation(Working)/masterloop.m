function [initial_power,reactivity,init1,init2,init3,init4,init5,init6,Toutold,TB,initial_iodine_pop,initial_xenon_pop,Tf] = masterloop(initial_power,reactivity,init1,init2,init3,init4,init5,init6,Toutold,TB,initial_iodine_pop,initial_xenon_pop,u,wierd,o,p)
    [initial_power,reactivity,init1,init2,init3,init4,init5,init6,Toutold,TB,initial_iodine_pop,initial_xenon_pop,Tf]=PRKEforloop(initial_power,reactivity,init1,init2,init3,init4,init5,init6,Toutold,TB,initial_iodine_pop,initial_xenon_pop,u,wierd,o,p);
    array1=[initial_power,reactivity,init1,init2,init3,init4,init5,init6,Toutold,TB,initial_iodine_pop,initial_xenon_pop,Tf];
    initial_power=array1(1);
    reactivity=array1(2);
    init1=array1(3);
    init2=array1(4);
    init3=array1(5);
    init4=array1(6);
    init5=array1(7);
    init6=array1(8);
    Toutold=array1(9);
    TB=array1(10);
    initial_iodine_pop=array1(11);
    initial_xenon_pop=array1(12);
    Tf=array1(13);
end 