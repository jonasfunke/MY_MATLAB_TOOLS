function [  ] = export_particle_stack_to_img(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

   [fname, pname] = uigetfile('*.mat', 'Select particle stack file');
    data = load([pname fname]);
    
    WriteImagic(data.particles, [pname fname(1:end-4) '_particles'])
    disp(['Particles written to: ' pname fname(1:end-4) '_particles.img'])
end

