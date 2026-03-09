
with open('SerraLINE/src/SerraLINE.chpl', 'r') as f:
    lines = f.readlines()

# Shift 244-556 left by 2 spaces
# Line numbers in the request are 1-based. 244 is index 243.
for i in range(243, 556):
    if lines[i].startswith('  '):
        lines[i] = lines[i][2:]
    elif lines[i].startswith(' '):
        lines[i] = lines[i][1:]

# Specifically fix the block that was reported as incorrect (290-337)
# Note: after the first shift, these lines still have 2 extra spaces.
# Original lines 290-337 (index 289-336)
for i in range(289, 337):
    if lines[i].startswith('  '):
        lines[i] = lines[i][2:]
    elif lines[i].startswith(' '):
        lines[i] = lines[i][1:]

with open('SerraLINE/src/SerraLINE.chpl', 'w') as f:
    f.writelines(lines)
