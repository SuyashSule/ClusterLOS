function refl_angle = reflectedAngle(center, point, incidentAngle)
radius_vector = point - center;
normal_vector = radius_vector/norm(radius_vector);
incident_vector = [cosd(incidentAngle) sind(incidentAngle)];
reflected_vector = incident_vector - 2*dot(incident_vector,normal_vector).*normal_vector;
refl_angle = -atan2d(-reflected_vector(2),reflected_vector(1));
end