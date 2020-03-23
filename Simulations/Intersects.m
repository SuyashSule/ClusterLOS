function intersects = Intersects(linkEndPt, BL_old_loc, BL_new_loc)
m1 = linkEndPt(2)/linkEndPt(1);
sign1 = (BL_old_loc(2)-m1*BL_old_loc(1))*(BL_new_loc(2)-m1*BL_new_loc(1));
m2 = (BL_old_loc(2)-BL_new_loc(2))/(BL_old_loc(1)-BL_new_loc(1));
c2 = BL_old_loc(2) - m2*BL_old_loc(1);
sign2 = (linkEndPt(2) - m2*linkEndPt(1) - c2)*(-c2);

if sign1 < 0 && sign2 < 0
    intersects = 1;
else
    intersects = 0;
end
end