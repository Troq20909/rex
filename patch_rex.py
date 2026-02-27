import re

content = open('rex_transpiler.py').read()

# Check if already patched
if "'using'" in content and "'size'" in content and "'on'" in content:
    print("Already patched - keywords present")
else:
    # Find the closing of the keywords set and add parameter keywords before it
    old = "'burgers_1d', 'navier_stokes_2d_simple', 'turbulence_stats'\n        }"
    new = ("'burgers_1d', 'navier_stokes_2d_simple', 'turbulence_stats',\n"
           "            # Parameter keywords used inside statements\n"
           "            'on', 'time', 'flux', 'eigenstates', 'using', 'mismatches', 'pam',\n"
           "            'length', 'dt', 'steps', 'velocity', 'viscosity', 'hot', 'cold',\n"
           "            'temp', 'compartments', 'basis', 'at', 'radius', 'of', 'size',\n"
           "            'all_results', 'to', 'hardening', 'for', 'this', 'model',\n"
           "        }")
    if old in content:
        content = content.replace(old, new)
        open('rex_transpiler.py', 'w').write(content)
        print("Patched successfully -", len(content), "bytes")
    else:
        print("Pattern not found - checking what's there:")
        # Show the keywords section
        idx = content.find('turbulence_stats')
        print(repr(content[idx:idx+50]))
