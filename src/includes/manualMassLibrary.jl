module massLibrary

import ..MasslistFunctions

CO2 = MasslistFunctions.createCompound(C=1, O=2)
NH3 = MasslistFunctions.createCompound(N=1, H=3)
PYRIDINE = MasslistFunctions.createCompound(N=1,H=5,C=5)
ACETONE = MasslistFunctions.createCompound(C=3,H=6,O=1)
APINENE = MasslistFunctions.createCompound(C=10, H=16)
BCARY = MasslistFunctions.createCompound(C=15, H=24)
ISOPRENE = MasslistFunctions.createCompound(C=5, H=8)
HEXANONE =  MasslistFunctions.createCompound(C=6, H=12, O=1)
PINONALDEHYDE = MasslistFunctions.createCompound(C=10, H=16, O=2)
PINONALDEHYDEPAN = MasslistFunctions.createCompound(C=10, H=15, O=6, N=1)
PINONICACID = MasslistFunctions.createCompound(C=10, H=16, O=3) #C10H16O3
PINICACID = MasslistFunctions.createCompound(C=9, H=14, O=4) #C9H14O4
ACETIC = MasslistFunctions.createCompound(C=2, H=4, O=2)
ACETICFRAG = MasslistFunctions.createCompound(C=2, H=2, O=1)
NH3 = MasslistFunctions.createCompound(N=1, H=3)
DMA = MasslistFunctions.createCompound(C=2, H=7, N=1)
TMA = MasslistFunctions.createCompound(C=3, H=9, N=1)
TMB = MasslistFunctions.createCompound(C=9, H=12)

ACETONITRILE = MasslistFunctions.createCompound(C=2, H=3, N=1)

OrgNitratNO = MasslistFunctions.createCompound(C=10, H=15, O=5, N=1)
NORPINONALDEHYDE = MasslistFunctions.createCompound(C=9, H=14, O=2)
NORPINONALDEHYDEPAN = MasslistFunctions.createCompound(C=9, H=13, O=6, N=1)
H3O = MasslistFunctions.createCompound(H=2, O=1)
H3OH2O = MasslistFunctions.createCompound(H=4, O=2)
H3OH2OH2O = MasslistFunctions.createCompound(H=6, O=3)

########## Sulphur ###############
DMS = MasslistFunctions.createCompound(C=2, H=6, S=1)
DMSO = MasslistFunctions.createCompound(C=2, H=6, S=1, O=1)
DMSO2 = MasslistFunctions.createCompound(C=2, H=6, S=1, O=2)
MSIA = MasslistFunctions.createCompound(C=1, H=4, S=1, O=2)
MSA = MasslistFunctions.createCompound(C=1, H=4, S=1, O=3)

########## BCARY Specific ########
#Fast
C15H22O2 = MasslistFunctions.createCompound(C=15, H=22, O=2)
C15H24O2 = MasslistFunctions.createCompound(C=15, H=24, O=2)
C15H24O3 = MasslistFunctions.createCompound(C=15, H=24, O=3)
C15H26O4 = MasslistFunctions.createCompound(C=15, H=26, O=4)
#Slow / Sticky
C14H20O3 = MasslistFunctions.createCompound(C=14, H=20, O=3)
C13H20O4 = MasslistFunctions.createCompound(C=13, H=20, O=4)
C14H22O4 = MasslistFunctions.createCompound(C=14, H=22, O=4)
C15H22O4 = MasslistFunctions.createCompound(C=15, H=22, O=4)
C15H24O4 = MasslistFunctions.createCompound(C=15, H=24, O=4)
C14H25NO4 = MasslistFunctions.createCompound(C=14, H=25, O=4, N=1)
C15H29NO6 = MasslistFunctions.createCompound(C=15, H=29, O=6, N=1)
#Nitrates
C15H23NO4 = MasslistFunctions.createCompound(C=15, H=23, O=4, N=1)
C13H19NO6 = MasslistFunctions.createCompound(C=13, H=19, O=6, N=1)
C15H25NO4 = MasslistFunctions.createCompound(C=15, H=25, O=4, N=1)
C15H21NO5 = MasslistFunctions.createCompound(C=15, H=21, O=5, N=1)
C15H23NO5 = MasslistFunctions.createCompound(C=15, H=23, O=5, N=1)
C15H23NO6 = MasslistFunctions.createCompound(C=15, H=23, O=6, N=1)
C15H17NO7 = MasslistFunctions.createCompound(C=15, H=17, O=7, N=1)

###### END BCARY ###############


###### NAPHTHA #################
NAPHTHA = MasslistFunctions.createCompound(C=10,H=8)
end
