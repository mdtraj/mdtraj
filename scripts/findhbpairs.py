
"""
Find the hydrogen bond pairs for a conformation
"""

acceptors = [('N', 'H'), ('NZ', 'HZ1'), ('NZ', 'HZ2'), ('NZ', 'HZ3'), ('OE1', 'HE1'), ('OE2', 'HE2'), ('OD1', 'HD1'), ('OD2', 'HD2'), ('ND1', 'HD1'), ('NE2', 'HE2'), ('NE2', 'HE21'), ('NE2', 'HE22'), ('ND2', 'HD21'), ('ND2', 'HD22'), ('OG', 'HG'), ('OG1', 'HG1'), ('OH', 'HH'), ('NH1', 'HH11'), ('NH1', 'HH12'), ('NH2', 'HH21'), ('NH2', 'HH22'), ('NE1', 'HE1'), ('SG', 'HG'), ('NE', 'HE')]
accHs = [tup[1] for tup in acceptors]
print accHs
donors = ['O', 'OD1', 'OD2', 'OE1', 'OE2', 'OG1', 'OG', 'OXT', 'OH', 'ND1', 'NE2', 'NZ', 'NE1', 'SG', 'SD']

donor_list = []
accH_list = []
for i, a in enumerate(t['AtomNames']):
    if a in donors:
        donor_list.append(i)

    if a in accHs:
        accH_list.append(i)
        
acc_list = []
to_pop = []
for i in accH_list:
    resID = t['ResidueID'][i]
    
    acc_name = [tup[0] for tup in acceptors if (tup[1] == t['AtomNames'][i])][0]
#    print acc_name, t['AtomNames'][i]
#    if t['AtomNames'][i] != 'H':
#        print t['AtomNames'][i], acc_name
#        break
#    acc_name = acc_name[0]
    try:
#        print t['AtomNames'][np.where(t['ResidueID'] == resID)]
        k = np.where((t['ResidueID'] == resID) & (t['AtomNames'] == acc_name))[0][0]
        acc_list.append(k)
    except:
        to_pop.append(i)

print to_pop
accH_list1 = []
for k in accH_list:
    if not k in to_pop:
        accH_list1.append(k)

print len(accH_list1), len(accH_list)
donor_list = np.array(donor_list)
acc_list = np.array(acc_list)
accH_list1 = np.array(accH_list1)
