lk_list = [];
for i = LinkeTurbidity(90-33,180-117,:)
lk_list = [lk_list,i];
end
plot(1:12, lk_list)
xlabel("Month")
ylabel("Linke Turbidity")
title("Linke Turbidity over the year")
% lk_list = [38   38   38   39   40   41   42   42   40   39   38   38]
