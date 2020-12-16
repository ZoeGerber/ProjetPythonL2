#################################################################################   SUJET  DM3    ################################################################################

#GERBER Zoé 20183203
#TIJANI Lobna 20181449
#MARTIN Fanny 20181851
#LEMERCIER Camille 20180595

###Enoncé###:
"""
A partir d'une séquence nucléotidique, vous recherchez tous les orfs dans les 3 phases directes,ainsi que celles sur le brin reverse complémentaire.
Vous donnerez les coordonnées de tous les orfs trouvés ainsi que leur séquence, leur longueur et leur traduction.
Commence avec ATG, se termine par TAG ou TAA ou TGA et doit tenir compte de la séquence à partir du premier caractère, puis la deuxième, puis de troisième
"""

#Définition d'un ORF sur wikipédia:
"""
En génétique moléculaire, un cadre de lecture ouvert, ou phase ouverte de lecture (open reading frame ou ORF en anglais)
est une partie d'un cadre de lecture susceptible d'être traduit en protéine ou en peptide.
C'est une suite de codons comprenant le codon start et un codon stop
Un codon d'initiation AUG du cadre de lecture ouvert — codon qui n'est pas nécessairement le 1er de celui-ci — peut indiquer le début de la traduction.
Le site de terminaison de la transcription est situé après le cadre de lecture ouvert, au-delà du codon stop.
La notion de cadre de lecture ouvert est très utilisée pour la prédiction de gènes.
Chaque séquence d'ADN peut contenir trois cadres de lecture décalés d'un nucléotide les uns par rapport aux autres (+1 ou -1).
Sur l'ADN, il peut y avoir transcription en ARN de l'un ou l'autre des 2 brins, ce qui conduit à un total de 6 cadres de lecture.
"""


##############Début du programme##############

#On rentre tout d'abord une séquence nucleotidique(ADN);

seqADN = input("Entrez une séquence nucléotidique : ")
print("\n")

#On vérifie si la séquence donnée est valide :

seqADN = seqADN.upper()
def verif (adn):
    i=0
    erreur = 0
    for i in adn :
        if not (( i == 'A' ) or ( i == 'T' ) or ( i == 'C' ) or ( i == 'G' )):
            erreur = erreur + 1
    return (erreur)

erreur = verif(seqADN)
if (verif(seqADN) == 0):#il n'y a donc aucune erreur dans la sequence donnée.
     print (seqADN, "est valide \n")
else:
    print (seqADN , "n'est pas une sequence valide")
    print ("il y a ", erreur , "erreur(s) dans cette séquence \n")

    
#Creation de la séquence de l'ADN complémentaire : 

def seqADNc (adn):
    i=0
    seqADNcompl = ""
    for i in adn :
        if (i == "A"):
            seqADNcompl = seqADNcompl + "T"
        elif (i == "T" ):
            seqADNcompl = seqADNcompl + "A"
        elif (i == "C"):
            seqADNcompl = seqADNcompl + "G"
        else :
            seqADNcompl = seqADNcompl + "C"
    return (seqADNcompl)

#Séquence du brin complémentaire

seqADNcompl=seqADNc(seqADN)
if (verif(seqADN) == 0):#affiche seulement si la séquence d'ADN donnée est valide.
    print ("la sequence d'ADN complémentaire à  notre sequence de base est : ", seqADNcompl,"\n")


#Création de la séquence reverse complémentaire : 

def reverseC(adn):
    i=0
    seqADNrc = ""
    for i in adn :
        if (i == "A"):
            seqADNrc = "A" + seqADNrc
        elif (i == "T" ):
            seqADNrc = "T" + seqADNrc
        elif (i == "C"):
            seqADNrc  = "C" + seqADNrc
        else :
            seqADNrc  = "G" + seqADNrc
    return (seqADNrc)   

#séquence du brin reverse complémentaire:
seqADNrc = reverseC(seqADNcompl)    
if (verif(seqADN) == 0):#affiche seulement si la séquence d'ADN donnée est valide.
    print ("la sequence d'ADN reverse complémentaire est : ", seqADNrc,"\n")



#Transcription de la sequence d'ADN en ARNmessager;

"""changer les T en U"""

def transcription (adn):
    i=0
    seqARN = ""
    for i in adn :
        if (i == "A"):
            seqARN = seqARN + "A"
        elif (i == "T" ):
            seqARN = seqARN + "U"
        elif (i == "C"):
            seqARN = seqARN + "C"
        else :
            seqARN = seqARN + "G"
    return (seqARN)

#séquence ARN du brin direct:

seqARN = transcription(seqADN)        
if (verif(seqADN) == 0):#affiche seulement si la séquence d'ADN donnée au debut est valide.
    print ("la sequence d'ARN transcrite est : ", seqARN ,"\n")


#sequence ARN du brin reverse complémentaire:

SeqARNrc = transcription(seqADNrc)
if (verif(seqADN) == 0):#affiche seulement si la séquence d'ADN donnée au debut est valide.
    print("la sequence ARN a partir du reverse complementaire:",SeqARNrc,"\n \n\n")

#Création du code génétique
#afin de permettre la traduction de l'ARNm

codeGene = {    "UUU" : "F",
                "UCU" : "S",
                "UAU" :	"Y",
                "UGU" :	"C",
                "UUC" :	"F",
                "UCC" :	"S",
                "UAC" :	"Y",
                "UGC" : "C",
                "UUA" :	"L",
                "UCA" : "S",
                "UAA" : "stop",
                "UGA":	"stop",
                "UUG":	"L",
                "UCG":	"S",
                "UAG":	"stop",
                "UGG":	"W",
                "CUU":	"L",
                "CCU":	"P",
                "CAU":	"H",
                "CGU":	"R",
                "CUC":	"L",
                "CCC":	"P",
                "CAC":	"H",
                "CGC":	"R",
                "CUA":	"L",
                "CCA":	"P",
                "CAA":	"Q",
                "CGA":	"R",
                "CUG":	"L",
                "CCG":	"P",
                "CAG":	"Q",
                "CGG":	"R",
                "AUU":	"I",
                "ACU":	"T",
                "AAU":	"N",
                "AGU":	"S",
                "AUC":	"I",
                "ACC":	"T",
                "AAC":	"N",
                "AGC":	"S",
                "AUA":	"I",
                "ACA":	"T",
                "AAA":	"K",
                "AGA":	"R",
                "AUG":	"M",
                "ACG":	"T",
                "AAG":	"K",
                "AGG":	"R",
                "GUU":	"V",
                "GCU":	"A",
                "GAU":	"D",
                "GGU":	"G",
                "GUC":	"V",
                "GCC":	"A",
                "GAC":	"D",
                "GGC":	"G",
                "GUA":	"V",
                "GCA":	"A",
                "GAA":	"E",
                "GGA":	"G",
                "GUG":	"V",
                "GCG":	"A",
                "GAG":	"E",
                "GGG":	"G",
           }

#######################################Traduction en proteine#################################################

taille=3
#création des fonctions qui permettent de commencer a traduire l'ARNm au premier codon start pour les trois différentes phases : 

    #Phase1#

def init1(arn):
    i=0
    seqARNtrad1=""
    for i in range(0,len(arn),3):
        codon=arn[i:i+3]
        if (codon=="AUG"):
            seqARNtrad1= arn[i:]
            break
    return(seqARNtrad1)


    #Phase 2#
    
def init2(arn):
    i=0
    seqARNtrad2=""
    for i in range(1,len(arn),3):
        codon=arn[i:i+3]
        if (codon=="AUG"):
            seqARNtrad2= arn[i:]
            break
    return (seqARNtrad2)


    #Phase 3#

def init3(arn):
    i=0
    seqARNtrad3=""
    for i in range(2,len(arn),3):
        codon=arn[i:i+3]
        if (codon=="AUG"):
            seqARNtrad3= arn[i:]
            break
    return (seqARNtrad3)


#création fonctions fin ORF, dans la phase 1(a partir de la sequence qui commence par un ATG):

def FinSeq(arn):
    i=0
    SeqStop1=""
    for i in range(0,len(arn),3):
        codon=arn[i:i+3]
        if ((codon=="UAG")or (codon=="UGA")or(codon=="UAA")):
            SeqStop1=arn[:i]
            break
        elif ((codon!="UAG")or (codon!="UGA")or(codon!="UAA")):
            SeqStop1=arn
    return (SeqStop1)


#Création de la fonction qui permet la traduction des séquences qui commence par un AUG : 

def tradC(arn):
    i=0
    prot = ""
    for i in range(0,len(arn)-taille+1,3) :
        codon = arn[i:i+taille]
        prot = prot + codeGene[codon]
        if ("stop" in prot):
            break
    return (prot)


#Création des 3 fonctions pour determiner les debuts des ORF:
"""permet de resortir la coordonnée du A du codon start"""

    #Phase 1

def debutORF1(arn):
    i=0
    debut1=0
    for i in range(0,len(arn),3):
        codon=arn[i:i+3]
        if (codon=="AUG"):
            debut1= i+1#car le i correspond a un chiffre en dessous de celui voulu.
            break
        elif(codon!="AUG"):
            debut1=0
    return (debut1)


    #Phase 2


def debutORF2(arn):
    i=0
    debut2=0
    for i in range(1,len(arn),3):
        codon=arn[i:i+3]
        if (codon=="AUG"):
            debut2= i+1#car le i correspond a un chiffre en dessous de celui voulu.
            break
        elif(codon!="AUG"):
            debut2=0
    return (debut2)

    #Phase3

def debutORF3(arn):
    i=0
    debut3=0
    for i in range(2,len(arn),3):
        codon=arn[i:i+3]
        if (codon=="AUG"):
            debut3= i+1#car le i correspond a un chiffre en dessous de celui voulu.
            break
        elif(codon!="AUG"):
            debut3=0
    return (debut3)

    
####fin des fonctions.

    ###Les proteines sur l'ARNm
    
#appel des fonctions pour sequences a traduire
""" création des sequences qui commence par un AUG dans les trois phases du brin direct"""

seqARN1 = init1(seqARN)#phase1
seqARN2 = init2(seqARN)#phase2
seqARN3 = init3(seqARN)#phase3

""" création des séquences  a traduire, qui commence a partir du codon start et qui finisse avant le codon stop dans les trois phases """

seqARNtrad1=FinSeq(seqARN1)#phase1
seqARNtrad2=FinSeq(seqARN2)#phase2
seqARNtrad3=FinSeq(seqARN3)#phase3

#appel des fonctions pour debut ORF

"""permet de ressortir les coordonnées du A du codon start """

debut1=0
debut1=debutORF1(seqARN)#phase1
debut2=0
debut2=debutORF2(seqARN)#phase2
debut3=0
debut3=debutORF3(seqARN)#phase3

#appel des fonctions pour fin ORF:
""" permet de recupérer les coordonnées du dernier nucleotides avant le codon stop trouvé dans les trois phases différentes
calcul de la taille de la sequence a traduire, on y ajoutera le nombre i, coordonnée du A afin de rajouter le debut de la séquence que nous avons precedement enlevé."""

fin1=0
fin1=len(seqARNtrad1)#phase1
fin2=0
fin2=len(seqARNtrad2)#phase2
fin3=0
fin3=len(seqARNtrad3)#phase3

### ARNm ###

if (verif(seqADN) == 0):#affiche seulement si la sequence d'ADN donnée est valide.
    print("-->les proteines sur l'ARNm\n \n")

    #Proteine1

    prot1 = tradC(seqARNtrad1)

    print("-Proteine en phase 1:\n")
    if len(prot1)==0:
        print("il n'y a pas de codon start dans la phase 1 du brin direct donc, pas de proteine dans cette phase.")

    elif len(prot1)!=0:
        print("sequence à traduire en phase 1 :",(seqARNtrad1))
        print ("La protéine 1 traduite selon le 1er ORF est : " , prot1)
        print ("La longueur de la proteine 1 est de ",len(prot1), "acide(s) aminé(s)")
        if ((fin1+debut1-1)!=len(seqARN)):
            print(" le debut de cet ORF est au ",debut1,"ème nucléotide, et la fin au ",(fin1+debut1-1), "ème nucléotide")
        else:
            print(" le debut de cet ORF est au ",debut1,"ème nucléotide, et la fin au ",len(seqADN), "ème nucléotide, (la fin de la sequence, il n'y a pas de codon stop)")

    print("\n")


#Proteine2

    prot2 = tradC(seqARNtrad2)
    print("-Proteine en phase 2:\n")

    if len(prot2)==0:
        print("il n'y a pas de codon start dans la phase 2 du brin direct donc, pas de proteine dans cette phase.")

    elif len(prot2)!=0:
        print("sequence à traduire en phase 2 :",(seqARNtrad2))
        print ("La protéine 2 traduite selon le 2eme ORF est : " , prot2)
        print ("La longueur de la proteine 2 est de ",len(prot2), "acide(s) aminé(s)")
        if ((fin2+debut2-1)!=len(seqARN)):
            print(" le debut de cet ORF est au ",debut2,"ème nucléotide, et la fin au ",(fin2+debut2-1), "ème nucléotide")
        else:
            print(" le debut de cet ORF est au ",debut2,"ème nucléotide, et la fin au ",len(seqADN), "ème nucléotide, (la fin de la sequence, il n'y a pas de codon stop)")

    print("\n")


#Proteine3

    prot3 = tradC(seqARNtrad3)
    print("-Proteine en phase 3:\n")

    if len(prot3)==0:
        print("il n'y a pas de codon start dans la phase 3 du brin direct donc, pas de proteine dans cette phase.")

    elif len(prot3)!=0:
        print("sequence à traduire en phase 3 :",(seqARNtrad3))
        print ("La protéine 3 traduite selon le 3eme ORF est : " , prot3)
        print ("La longueur de la proteine 3 est de ",len(prot3), "acide(s) aminé(s)")
        if ((fin3+debut3-1)!=len(seqARN)):
            print(" le debut de cet ORF est au ",debut3,"ème nucléotide, et la fin au ",(fin3+debut3-1), "ème nucléotide")
        else:
            print(" le debut de cet ORF est au ",debut3,"ème nucléotide, et la fin au ",len(seqADN), "ème nucléotide, (la fin de la sequence, il n'y a pas de codon stop)")

    print("\n\n\n")

 ##Les proteines sur l'ADN reverse complementaire##

#Transcription du reverse complementaire#


#Appel des fonctions pour creer la sequence a traduire
""" commencer a un AUG dans les trois différentes phases"""

seqARN1rc = init1(SeqARNrc)#phase1
seqARN2rc = init2(SeqARNrc)#phase2
seqARN3rc = init3(SeqARNrc)#phase3

""" finir avant le codon stop pour les trois phases également"""

seqARNtrad1rc=FinSeq(seqARN1rc)#phase1
seqARNtrad2rc=FinSeq(seqARN2rc)#phase2
seqARNtrad3rc=FinSeq(seqARN3rc)#phase3

#appel des fonctions pour debut ORF

debut1rc=0
debut1rc=debutORF1(SeqARNrc)#phase1
debut2rc=0
debut2rc=debutORF2(SeqARNrc)#phase2
debut3rc=0
debut3rc=debutORF3(SeqARNrc)#phase3

#appel des fonctions pour fin ORF:


fin1rc=0
fin1rc=len(seqARNtrad1rc)#phase1


fin2rc=0
fin2rc=len(seqARNtrad2rc)#phase2

fin3rc=0
fin3rc=len(seqARNtrad3rc)#phase3


    ## ARN reverse complémentaire ##

if (verif(seqADN) == 0):#si la sequence d'ADN donnée n'est pas valide, cela n'affichera rien.
    print("-->les proteines sur l'ARNm du reverse complémentaire:\n \n")

#Proteine1
    
    print("-Proteine en phase 1:\n")

    prot1rc = tradC(seqARNtrad1rc)
    
    if len(prot1rc)==0:
        print("il n'y a pas de codon start dans la phase 1 du brin reverse complémentaire donc, pas de proteine dans cette phase.")

    elif len(prot1rc)!=0:
        print("sequence à traduire en phase 1 :",(seqARNtrad1rc))
        print ("La protéine 1 traduite selon le 1er ORF est : " , prot1rc)
        print ("La longueur de la proteine 1 est de ",len(prot1rc), "acide(s) aminé(s)")
        if ((fin1rc+debut1rc-1)!=len(seqARN)):
            print(" le debut de cet ORF est au ",debut1rc,"ème nucléotide, et la fin au ",(fin1rc+debut1rc-1), "ème nucléotide")
        else:
            print(" le debut de cet ORF est au ",debut1rc,"ème nucléotide, et la fin au ",len(seqADN), "ème nucléotide, (la fin de la sequence, il n'y a pas de codon stop)")

    print("\n")

#Proteine2
    
    
    print("-Proteine en phase 2:\n")
    prot2rc = tradC(seqARNtrad2rc)
  
    if len(prot2rc)==0:
        print("il n'y a pas de codon start dans la phase 2 du brin reverse complémentaire donc, pas de proteine dans cette phase.")

    elif len(prot2rc)!=0:
        print("sequence à traduire en phase 2 :",(seqARNtrad2rc))
        print ("La protéine 2 traduite selon le 2eme ORF est : " , prot2rc)
        print ("La longueur de la proteine 2 est de ",len(prot2rc), "acide(s) aminé(s)")
        if ((fin2rc+debut2rc-1)!=len(seqARN)):
            print(" le debut de cet ORF est au ",debut2rc,"ème nucléotide, et la fin au ",(fin2rc+debut2rc-1), "ème nucléotide")
        else:
            print(" le debut de cet ORF est au ",debut2rc,"ème nucléotide, et la fin au ",len(seqADN), "ème nucléotide, (la fin de la sequence, il n'y a pas de codon stop)")


    print("\n")


#Proteine3

    print("-Proteine en phase 3:\n")
    prot3rc = tradC(seqARNtrad3rc)
    


    if len(prot3rc)==0:
        print("il n'y a pas de codon start dans la phase 3 du brin reverse complémentaire donc, pas de proteine dans cette phase.")

    elif len(prot3rc)!=0:
        print("sequence à traduire en phase 3 :",(seqARNtrad3rc))
        print ("La protéine 3 traduite selon le 3eme ORF est : " , prot3rc)
        print ("La longueur de la proteine 3 est de ",len(prot3rc), "acide(s) aminé(s)")
        if ((fin3rc+debut3rc-1)!=len(seqARN)):
            print(" le debut de cet ORF est au ",debut3rc,"ème nucléotide, et la fin au ",(fin3rc+debut3rc-1), "ème nucléotide")
        else:
            print(" le debut de cet ORF est au ",debut3rc,"ème nucléotide, et la fin au ",len(seqADN), "ème nucléotide, (la fin de la sequence, il n'y a pas de codon stop)")
    print("\n\n")

