%function parsave(fname, x, y, z, w, u, v, p, q)
function parsave(fname, x, y, z, w, u)

% save files in parfor loop to avoid workspace transparency error
% fname is the filename (as a string)
% x, y, z, w, u, v, p, q are workspace variables

  %renamevars2save(x,y,z,w,u,v,p,q,inputname(2),inputname(3),inputname(4),inputname(5), inputname(6), inputname(7),inputname(8),inputname(9))
  %save(fname, inputname(2), inputname(3), inputname(4), inputname(5), inputname(6), inputname(7),inputname(8),inputname(9))
  
  renamevars2save(x,y,z,w,u,inputname(2),inputname(3),inputname(4),inputname(5),inputname(6))
  save(fname, inputname(2), inputname(3), inputname(4), inputname(5),inputname(6))
end

%function renamevars2save(x,y,z,w,u,v,p,q,x_name,y_name,z_name,w_name,u_name,v_name, p_name, q_name)
function renamevars2save(x,y,z,w,u,x_name,y_name,z_name,w_name,u_name)
    assignin('caller', x_name, x)
    assignin('caller', y_name, y)
    assignin('caller', z_name, z)
    assignin('caller', w_name, w)
    assignin('caller', u_name, u)
%     assignin('caller', v_name, v)
%     assignin('caller', p_name, p)
%     assignin('caller', q_name, q)    
end