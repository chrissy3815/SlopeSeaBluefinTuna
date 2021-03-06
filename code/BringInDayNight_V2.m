
addpath /Users/chrissy/JointProgram/SlopeSea/2016_data/SlopeSeaBluefinTuna/data
addpath /Users/chrissy/JointProgram/SlopeSea/2016_data/SlopeSeaBluefinTuna/code

%kk=find(b(k)<rs(1) & b(k)>rs(2))
T=readtable('GU1608HB1603Event.csv');
T.TimeDecimal=hour(T.EVENT_DATE)+minute(T.EVENT_DATE)/60;
T.DayNight(:)="Unk";
T.Sunrise(:)=nan;
T.Sunset(:)=nan;
for n=1:height(T)
    [rs,t,d,z,a,r]=suncycle(T.LATITUDE(n),T.LONGITUDE(n),T.EVENT_DATE(n));
    T.Sunrise(n)=rs(1);
    T.Sunset(n)=rs(2);
    if rs(2)<rs(1)
        if T.TimeDecimal(n)>rs(2) & T.TimeDecimal(n)<rs(1)
            T.DayNight(n)="Night";
        end
        if T.TimeDecimal(n)<rs(2) | T.TimeDecimal(n)>rs(1)
            T.DayNight(n)="Day";
        end
    end
    if rs(1)<rs(2)
        if T.TimeDecimal(n)>rs(1) & T.TimeDecimal(n)<rs(2)
            T.DayNight(n)="Day";
        end
        if T.TimeDecimal(n)<rs(1) | T.TimeDecimal(n)>rs(2)
            T.DayNight(n)="Night";
        end
    end
    
end

% save as a csv
writetable(T,'/Users/chrissy/JointProgram/SlopeSea/2016_data/SlopeSeaBluefinTuna/data/GU1608HB1603Event_withDayNight.csv')

