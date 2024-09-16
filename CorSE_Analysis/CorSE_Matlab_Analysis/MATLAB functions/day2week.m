function week = day2week(day)

a = 1;
WeekVec = [0 0 0 0 0 0 0];

for i = 1:20
    
    WeekVec = [WeekVec a a a a a a a];
    
    a = a + 1;
    
end

week = WeekVec(day);

