function quantitiesForExpansion(sc,rp)
%GETQUANTITIES get all quantities  that can be provided by this interface
%   This method is called by the class RieszProjection, if it is
%   constructed with an instance of the class 'Scattering'.

sc.parent = rp; 
quantities = struct; 
qNames = fieldnames(sc.quantities);
for nm = qNames.'
    if ~sc.quantities.(nm{1}).logical(3), continue; end
    quantities.(nm{1}).quadratic = sc.quantities.(nm{1}).logical(1);
    quantities.(nm{1}).hidden = sc.quantities.(nm{1}).logical(2);
    quantities.(nm{1}).m = sc.quantities.(nm{1}).m;
    quantities.(nm{1}).parent = sc.quantities.(nm{1}).parent;
    if isfield(sc.quantities.(nm{1}),'plot')
        quantities.(nm{1}).plot = sc.quantities.(nm{1}).plot;
    end
    switch nm{1}
        case 'DipoleEmission'
            quantities.DipoleEmission.evaluate = @dipoleEmission;
            ndk.parents = {'DipoleEmission'};
            ndk.evaluate = @(dpe,w0)dpe./dipoleBulkEmission(sc,w0);
            ndk.quadratic = false;
            quantities.NormalizedDecayRate = ndk;
        case 'FarFieldIntegral'
            pce.quadratic = true;
            pce.parents = {'FarFieldIntegral'};
            pce.evaluate = @(x)x(:,1,:)./sum(x(:,2,:),3);
            quantities.PhotonCollectionEfficiency = pce;
            if sc.pointSource
                dpce.quadratic = true;
                dpce.parents = {'FarFieldIntegral' 'DipoleEmission'};
                dpce.evaluate = @(rad,dpe)rad(:,1,:)./real(sum(dpe,3));
                quantities.DipolePowerCollectionEfficiency = dpce;
            end
            rad.quadratic = true;
            rad.parents = {'FarFieldIntegral'};
            rad.evaluate = @(x)x(:,1,:);
            quantities.Radiation = rad;
    end
end
rp.quantities = quantities;

    function dpe = dipoleEmission(dpe,w0)
        if ~isreal(dpe)
            dpe(:,:,end) = dpe(:,:,end) + dipoleBulkEmission(sc,w0);
        end
    end
end