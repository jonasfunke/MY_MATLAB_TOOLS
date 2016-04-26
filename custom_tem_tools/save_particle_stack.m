function [] = save_particle_stack(data, fileloc)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    angles = data.angles;
    history = data.history;
    name = data.name;
    particles = data.particles;
    save(fileloc, 'angles', 'history', 'name', 'particles')
    disp(['Stack written to: ' fileloc])
end

