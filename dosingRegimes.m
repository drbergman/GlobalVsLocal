function EVENTS = dosingRegimes(censor_date,DoW_start,...
    aFGFR3_start_simday,aFGFR3_days_between,...
    n_doses_aFGFR3,break_after,break_length)

n_subcohorts = length(aFGFR3_days_between);
EVENTS = cell(n_subcohorts,1);


for i = 1:n_subcohorts
    
    event_schedule = [aFGFR3_start_simday(i)+aFGFR3_days_between(i)*(0:n_doses_aFGFR3(i)-1)';... % aFGFR3 doses
        censor_date]; % censoring
    
    event_ind = [2*ones(n_doses_aFGFR3(i),1);...
        Inf];
    
    EVENTS{i} = sortrows([event_schedule,event_ind]);
    EVENTS{i} = avoidWeekends(EVENTS{i},DoW_start,censor_date);
    
    if ~isempty(break_after)
        for drug_ind = 1:2
            first_drug = EVENTS{i}(find(EVENTS{i}(:,2)==drug_ind,1),1);
            later_dosings = (EVENTS{i}(:,1)>=(first_drug+break_after(i))) & (EVENTS{i}(:,2)==drug_ind);
            EVENTS{i}(later_dosings,1) = EVENTS{i}(later_dosings,1) + break_length(i);
        end
        EVENTS{i} = avoidWeekends(EVENTS{i},DoW_start,censor_date);
    end
end