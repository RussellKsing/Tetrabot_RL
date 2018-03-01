function [f,g]=get_grad(x,getf,getg)
    f=getf(x);
    g=getg(x);
end