function EVENTS = avoidWeekends(EVENTS,start_day,censor_date)

% make all dosings occur on weekdays

weekend_dose = find(any(mod(start_day + EVENTS(:,1),7) == [5,6],2) & any(EVENTS(:,2)==[1,2],2),1); % earliest therapy scheduled for a weekend

while ~isempty(weekend_dose)
    future_sametype_doses = (EVENTS(:,2)==EVENTS(weekend_dose,2)) & (1:size(EVENTS,1))'>=weekend_dose; % all therapies of the same type as the next weekend therapy and occurring on or after that date; these will all be pushed back by 1 or 2 days depending on whether the therapy is scheduled for Sunday or Saturday, resectively
    
    EVENTS(future_sametype_doses,1) = EVENTS(future_sametype_doses,1) + mod(-start_day-EVENTS(weekend_dose,1),7);
    
    weekend_dose = find(any(mod(start_day + EVENTS(:,1),7) == [5,6],2) & any(EVENTS(:,2)==[1,2],2),1); % next remaining therapy scheduled for a weekend
end

EVENTS = sortrows(EVENTS);

EVENTS(EVENTS(:,1)>=censor_date & EVENTS(:,2)~=Inf,:) = [];