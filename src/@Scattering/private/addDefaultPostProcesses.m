function  addDefaultPostProcesses(sc)
%ADDDEFAULTQUANTITIES adds default quantities for expansions

reference = {'referenceDipoleEmission'};
sc.addPostProcess('DipoleEmission',sc.keys,reference);
sc.addPostProcess('PointEvaluation',false)

keys = struct('interfaceType','ExteriorDomain');
keys_ref = keys; 
keys_ref.outputQuantity = 'ElectromagneticFieldEnergyFlux';
reference = {'referenceFluxIntegration' keys_ref};
sc.addPostProcess('ElectromagneticFieldEnergyFlux',keys,reference,true)

keys_ref = struct('outputQuantity','ElectromagneticFieldAbsorption');
reference = {'referenceDensityIntegration' keys_ref};
sc.addPostProcess('ElectromagneticFieldAbsorption',reference,true);

keys_ref = struct('outputQuantity','ElectricFieldEnergy');
reference = {'referenceDensityIntegration' keys_ref};
sc.addPostProcess('ElectricFieldEnergy',reference,true);

sc.addPostProcess('FarFieldIntegral','FarField',true)
sc.addPostProcess('RadiationPattern','FarField',true)

sc.addPostProcess('FourierTransform',false);
sc.addPostProcess('ModeCoupling',true);

sc.addPostProcess('CartesianExport',false);
sc.addPostProcess('Superposition',false);
end

