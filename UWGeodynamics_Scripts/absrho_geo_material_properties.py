# core UW bit
import UWGeodynamics as GEO
u = GEO.UnitRegistry

import geo_model_properties as modprop

alpha   = 3.0e-5 / u.degK

# Density is defined using the dimensionless relative density
um_surf_density  = 3400 * u.kilogram / u.metre**3 # at surface
abs_density      = um_surf_density * (1 - alpha * (modprop.Tint - modprop.Tsurf))
ref_density      = um_surf_density * (modprop.Tint-modprop.Tsurf) * alpha

eclogite_surf_density = 3500 * u.kilogram / u.meter**3
eclogite_density = eclogite_surf_density * ( 1 - alpha * (modprop.Tint - modprop.Tsurf))
eclogite_density2=eclogite_surf_density * ( 1 - alpha * (modprop.Tint - modprop.Tsurf))
abs_density, ref_density, eclogite_density

modprop.Tint

modprop.Tsurf

# Material properties
um = {
    'name'     : 'upper mantle', # really the 'upper mantle'
    'index'    : 0,
    'viscosity': 1. * u.pascal * u.second * 1e20,
    'density'  : abs_density
}

lm = {
    'name'     : 'lower mantle',
    'index'    : 1,
    'viscosity': 30. * u.pascal * u.second * 1e20,
    'density'  : 3407.9742* u.kg / u.m**3
}


# 80Ma oceanic lithosphere
subplate1 = {
    'name'     : 'oceanic plate 1',
    'index'    : 3,
    'viscosity': 1.000e+25 * u.pascal * u.second,
    'density'  : 3242.49 * u.kg / u.m**3,
    'cohesion' : 1.250e+01 * u.megapascal,
    'cohesion2': 6.250e+00 * u.megapascal,
}

subplate1_phase = {
    'name'     : 'oceanic plate 1 after phase change',
    'index'    : 2,
    'viscosity': 1.000e+25 * u.pascal * u.second,
    'density'  : eclogite_density,
    'cohesion' : 1.250e+01 * u.megapascal,
    'cohesion2': 6.250e+00 * u.megapascal,
}

subplate2 = {
    'name'     : 'oceanic plate 2',
    'index'    : 4,
    'viscosity':  1.000e+25 * u.pascal * u.second,
    'density'  :  3347.5 * u.kg / u.m**3,
    'cohesion' :  66.94385 * u.megapascal,
    'cohesion2':  33.4719 * u.megapascal,
}


subplate3 = {
    'name'     : 'oceanic plate 3',
    'index'    : 5,
    'viscosity': 1.929e+24 * u.pascal * u.sec,
    'density'  : 3318.5 * u.kg / u.m**3,
    'cohesion' : 120.3433 * u.megapascal,
    'cohesion2': 60.1716 * u.megapascal,
}

subplate4 = {
    'name'     : 'oceanic plate 4',
    'index'    : 6,
    'viscosity': 9.64083e+21 * u.pascal * u.second,
    'density'  : 3297.0 * u.kg / u.m**3,
    # no yielding
}

# +
# cratonic lithosphere
craton1 = {
    'name'     : 'craton1',
    'index'    : 11,
    'viscosity': 5.00000e+23 * u.pascal * u.second,
    'density'  : 2.98662e+03  * u.kg / u.m**3,
    'cohesion' : 1.3000e+02 * u.megapascal,
    'cohesion2': 6.500e+01 * u.megapascal,
}

craton2 = {
    'name'     : 'craton2',
    'index'    : 12,
    'viscosity': 1.67484e+23 * u.pascal * u.second,
    'density'  : 3.30100e+03 * u.kilogram / u.m**3,
    'cohesion' : 2.91998e+02 * u.megapascal,
    'cohesion2': 1.45999e+02 * u.megapascal,
}
# -

# transitional lithosphere
trans1 = {
    'name'     : 'trans1',
    'index'    : 9,
    'viscosity': 5.00000e+23 * u.pascal * u.second,
    'density'  : 3.00529e+03 * u.kg / u.m**3,
    'cohesion' : 4.000e+01 * u.megapascal,
    'cohesion2': 2.000e+01 * u.megapascal,
}
trans2 = {
    'name'     : 'trans2',
    'index'    : 10,
    'viscosity': 1.50142e+23 * u.pascal * u.second,
    'density'  : 3.30122e+03 * u.kg / u.m**3,
    'cohesion' : 1.500e+02 * u.megapascal,
    'cohesion2': 7.500e+01 * u.megapascal,
}

# weak back arc material properties
backArc1 = {
    'name'     : 'backArc1',
    'index'    : 7,
    'viscosity': 4.97782e+23 * u.pascal * u.second,
    'density'  : 3.10862e+03 * u.kg/u.m**3,
    'cohesion' : 1.250e+01 * u.megapascal,
    'cohesion2': 6.250e+00 * u.megapascal,
}
backArc2 = {
    'name'     : 'backArc2',
    'index'    : 8,
    'viscosity': 1.71934e+22 * u.pascal * u.second,
    'density'  : 3.28283e+03 * u.kg/u.m**3,
    'cohesion' : 2.500e+01 * u.megapascal,
    'cohesion2': 1.250e+01 * u.megapascal,
}

# +
# assume ribbon and buoyant strip have cratonic material properties
ribbon1 = {
    'name'     : 'ribbon_1',
    'index'    : 13,
    'viscosity': 1e25 * u.pascal * u.second,
    'density'  : 2840.5709  * u.kg / u.m**3,
    'cohesion' : 12.5 * u.megapascal,
    'cohesion2': 6.25 * u.megapascal,
}

ribbon2 = {
    'name'     : 'ribbon_2',
    'index'    : 14,
    'viscosity': 1e25 * u.pascal * u.second,
    'density'  : 3280.8644 * u.kg / u.m**3,
    'cohesion' : 57.944 * u.megapascal,
    'cohesion2': 6.25 * u.megapascal,
}
# -

buoyStrip = {
    'name'     : 'buoyStrip',
    'index'    : 15,
    'viscosity': 1e25 * u.pascal * u.second,    # strong
    'density'  : 2.98662e+03  * u.kg / u.m**3,         # assume cratonic density
}

# define material list
material_list = [ um, lm,
                 subplate1, subplate1_phase, subplate2, subplate3, subplate4,
                 backArc1, backArc2,
                 trans1, trans2,
                 craton1, craton2,
                 ribbon1,ribbon2,
                 buoyStrip ]
