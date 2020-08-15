import MMesh5
v,f=MMesh5.cylinder_with_spherical_caps(1,2.5,3,15)
e=MMesh5.get_edge(f)
#print(e)
#print(v)
print(len(v),len(f),len(e))
print(MMesh5.euler_number(v,f))