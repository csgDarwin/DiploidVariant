import sys
import os
import Bio


class CladeAlignDB:
    def __init__(self,DBase):
        self.db=DBase
        self.nameList=[]
        self.ref_list=[]
        self.numref_species=0
        #print("newAlignDB created")
        if self.db=="20birds":
            self.OrderList = ["oscTAEGUT","oscVIDCHA","oscCORHAW","oscCORMON","oscAGEPHO","oscCATUST","oscHIRRUS","parMELUND","parSTRHAB","capCALANA","subCHILAN",
                         "falFALNAU","falFALPER","falFALBIA","falFALCHE","falFALRUS","eagACCGEN","eagAQUCHR","capAPUAPU","outGALGAL"]
            self.NameDict={"oscTAEGUT":"Zebra_finch","oscVIDCHA":"Village_indigobird","oscCORHAW": "Hawaiian_crow", "oscCORMON":"New_Caledonian_crow",
                       "oscAGEPHO":"Red-winged_blackbird","oscCATUST":"Swainsons_thrush","oscHIRRUS":"Barn_swallow","subCHILAN":"Lanced-tailed_manakin",
                       "parMELUND":"Budgerigar","parSTRHAB":"Kakapo","capCALANA":"Annas_hummingbird","falFALNAU":"Lesser_kestrel","falFALPER":"Peregrine_falcon",
                       "falFALBIA":"Lanner_falcon","falFALCHE":"Saker_falcon","falFALRUS":"Gyrfalcon","eagACCGEN":"Northern_goshawk","eagAQUCHR":"European_golden_eagle",
                       "capAPUAPU":"Common_swift","outGALGAL":"Chicken"}
            self.num_species=20
            self.numref_species=10
            self.ref_list = ["oscTAEGUT","oscVIDCHA","oscCORHAW","oscCORMON","oscAGEPHO","oscCATUST","oscHIRRUS","parMELUND","parSTRHAB","capCALANA"]
            #self.ApeList = ["hg38","panTro5","panPan2","gorGor5","ponAbe2","nomLeu3"]
        elif self.db=="45way":
            self.OrderList = ["priHomSap38","priHomSapT2T","oscTAEGUT","oscVIDCHA","oscCORHAW","oscCORMON","oscAGEPHO","oscCATUST","oscHIRRUS","parMELUND",
                         "parSTRHAB","capCALANA",
                         "priPanTro","priPanPan","priGorGor","priCalJac","priAotNan","priCarSyr","rodMusMus","rodRatNor",
                         "rodArvNil","dolTurTru","dolPhoSin","whaBalMus","rumBosTau","rumOviAri","rumCerEla","batRhiFer","batPhyDis","mamDasNov","mamOrnAna",
                         "subCHILAN","falFALNAU","falFALPER","falFALBIA","falFALCHE","falFALRUS","eagACCGEN","eagAQUCHR","capAPUAPU","outGALGAL","turTraScr","lizPogVit",
                         "lizLacAgi","ampRanTem"]
            self.NameDict={"oscTAEGUT":"Zebra_finch","oscVIDCHA":"Village_indigobird","oscCORHAW": "Hawaiian_crow", "oscCORMON":"New_Caledonian_crow",
                       "oscAGEPHO":"Red-winged_blackbird","oscCATUST":"Swainsons_thrush","oscHIRRUS":"Barn_swallow","subCHILAN":"Lanced-tailed_manakin",
                       "parMELUND":"Budgerigar","parSTRHAB":"Kakapo","capCALANA":"Annas_hummingbird","falFALNAU":"Lesser_kestrel","falFALPER":"Peregrine_falcon",
                       "falFALBIA":"Lanner_falcon","falFALCHE":"Saker_falcon","falFALRUS":"Gyrfalcon","eagACCGEN":"Northern_goshawk","eagAQUCHR":"European_golden_eagle",
                       "capAPUAPU":"Common_swift","outGALGAL":"Chicken",
                       "priHomSapT2T":"Human_hs1","priHomSap38":"Human_hg38","priPanTro":"Chimpanzee","priPanPan":"Bonobo","priGorGor":"Gorilla","priCalJac":"Marmoset",
                       "priAotNan":"Ma's_Night_Monkey","priCarSyr":"Philippine_Tarsier","rodMusMus":"House_Mouse","rodRatNor":"Norway_Rat",
                         "rodArvNil":"Nile_Rat","dolTurTru":"Bottlenose_Dolphin","dolPhoSin":"Vaquita","whaBalMus":"Blue_Whale","rumBosTau":"Cow","rumOviAri":"Sheep",
                         "rumCerEla":"Red_Deer","batRhiFer":"Greater_horseshoe_bat","batPhyDis":"Pale_spear-nosed_bat","mamDasNov":"nine-banded_armadillo","mamOrnAna":"Platypus",
                         "subCHILAN":"lanced-tailed_manakin","falFALNAU":"lesser_kestrel","falFALPER":"peregrine_falcon","falFALBIA":"Lanner_falcon","falFALCHE":"saker_falcon",
                         "falFALRUS":"gyrfalcon","eagACCGEN":"Northerna_goshawk","eagAQUCHR":"European_golden_eagle","capAPUAPU":"common_swift","outGALGAL":"Domestic_chicken",
                         "turTraScr":"Red-eared_Slider","lizPogVit":"Central_bearded_dragon","lizLacAgi":"Sand_lizard","ampRanTem":"European_common_frog"
                       }
            self.num_species=45
            self.numref_species=12
            self.ref_list =["priHomSap38","priHomSapT2T","oscTAEGUT","oscVIDCHA","oscCORHAW","oscCORMON","oscAGEPHO","oscCATUST","oscHIRRUS","parMELUND","parSTRHAB","capCALANA"]
        elif self.db=="gene_wide":
            self.OrderList = ["pHomSap38","pHomSapT2T","pPanPan","pPanTro","pGorGor","pPonPyg","pPonAbe","pSymSyn","pCalJac","pNycCou","pLemCat","mCynVol",
                         "mOchPri","mSciCar","mJacJac","mPerMan","mOnyTor","mChiNiv","mArvAmp","mAcoRus","mRatNor","mArvNil","mApoSyl","mMusMus",
                         "mDelDel","mTurTru","mGloMel","mLagAlb","mOrcOrc","mPhoSin","mMesDen","mKogBre","mEubGla","mBalMus","mBalAcu","mHipAmp","mBosTau",
                         "mCerEla","mSusScr","mCamDro","mEquCab","mDicBic","mManPen","mLynCan","mNeoNeb","mCanLup","mZalCal","mMelMel","mLutLut","mMusErm",
                         "mRhiFer","mRouAeg","mPhyDis","mDesRot","mMolMol","mPipKuh","mMyoDau","mMyoMyo","mEriEur","mSunEtr","mSorAra","mLoxAfr","mEleMax",
                         "mDasNov","mChoDid","mMonDom","mDroGli","mSarHar","mTriVul","mPhaCin","mTacAcu","mOrnAna","bTaeGut","bVidCha","bCorHaw","bCorMon",
                         "bMelGeo","bMelMel","bAmmCau","bMolAte","bHaeMex","bSerCan","bCinCin","bCatUst","bPoeAtr","bHirRus","bChiLan","bMelUnd","bFalNau",
                         "bFalPer","bFalBia","bFalRus","bFalChe","bAquChr","bAccGen","bColStr","bPogPus","bDryPub","bRisTri","bChrRid","bGruAme","bGavSte",
                         "bCalAnn","bApuApu","bCucCan","bCygOlo","bAytFul","bAnaPla","bLagMut","bGalGal","bRhePen","rAllMis","rGopFla","rMalTer","rDerCor",
                         "rCheMyd","rCarCar","rEulEur","rThaEle","rCanAsp","rElgMul","rRhiFlo","rZooViv","rLacAgi","rPodRaf","aMicUni","aGeoSer","aRhiBiv",
                         "aBomBom","aPelFus","aRanTem","aBufBuf","aHylSar","aPseCor"]
            self.NameDict={
                        "pHomSap38":"hg38","pHomSapT2T":"hs1","pPanPan":"Bonobo","pPanTro":"Chimpanzee","pGorGor":"Gorilla","pPonPyg":"Bornean orangutan",
                        "pPonAbe":"Sumatran orangutan","pSymSyn":"Siamang gibbon","pCalJac":"Common marmoset","pNycCou":"Sunda slow loris","pLemCat":"Ring-tailed lemur",
                        "mCynVol":"Philippine flying lemur","mOchPri":"American pika","mSciCar":"Grey squirrel","mJacJac":"Lesser egyptian jerboa","mPerMan":"Deer mouse",
                        "mOnyTor":"Southern grasshopper mouse","mChiNiv":"European snow vole","mArvAmp":"European water vole","mAcoRus":"Golden spiny mouse",
                        "mRatNor":"Norway rat","mArvNil":"Nile rat","mApoSyl":"European woodmouse","mMusMus":"House mouse","mDelDel":"Saddleback dolphin",
                        "mTurTru":"Bottlenose dolphin","mGloMel":"Long-finned pilot whale","mLagAlb":"White-beaked dolphin","mOrcOrc":"Killer whale",
                        "mPhoSin":"Vaquita","mMesDen":"Blainville's beaked whale","mKogBre":"Pygmy sperm whale","mEubGla":"North atlantic right whale",
                        "mBalMus":"Blue whale","mBalAcu":"Common minke whale","mHipAmp":"Common hippopotamus","mBosTau":"Brahman cow","mCerEla":"Red deer",
                        "mSusScr":"Pig","mCamDro":"Dromedary camel","mEquCab":"Arabian horse","mDicBic":"Southern black rhinoceros","mManPen":"Chinese pangolin",
                        "mLynCan":"Canada lynx","mNeoNeb":"Clouded leopard","mCanLup":"German shepherd dog","mZalCal":"California sea lion","mMelMel":"European badger",
                        "mLutLut":"Eurasian otter","mMusErm":"Stout","mRhiFer":"Greater horseshoe bat","mRouAeg":"Egyptian fruit bat","mPhyDis":"Pale spear-nosed bat",
                        "mDesRot":"Common vampire bat","mMolMol":"Pallas's mastiff bat","mPipKuh":"Kuhl's pipistrelle","mMyoDau":"Daubenton's bat",
                        "mMyoMyo":"Greater mouse-eared bat","mEriEur":"European hedgehog","mSunEtr":"Etruscan shrew","mSorAra":"European shrew",
                        "mLoxAfr":"African elephant","mEleMax":"Asian elephant","mDasNov":"Nine-banded armadillo","mChoDid":"Linnaeus's two-toed sloth",
                        "mMonDom":"Gray short-tailed opossum","mDroGli":"Monito del montes","mSarHar":"Tasmanian devil","mTriVul":"Common brushtail possum",
                        "mPhaCin":"Koala","mTacAcu":"Short-beaked echidna","mOrnAna":"Platypus","bTaeGut":"Zebra finch","bVidCha":"Village indigobird",
                        "bCorHaw":"Hawaiian crow","bCorMon":"New caledonian crow","bMelGeo":"Swamp sparrow","bMelMel":"Song sparrow","bAmmCau":"Saltmarsh sparrow",
                        "bMolAte":"Brown-headed cowbird","bHaeMex":"House finch","bSerCan":"Canary","bCinCin":"White-throated dipper",
                        "bCatUst":"Swainson's thrush","bPoeAtr":"Black-capped chickadee","bHirRus":"Barn swallow","bChiLan":"Lanced-tailed manakin",
                        "bMelUnd":"Budgerigar","bFalNau":"Lesser kestrel","bFalPer":"Peregrine falcon","bFalBia":"Lanner falcon","bFalRus":"Gyrfalcon",
                        "bFalChe":"Saker falcon","bAquChr":"European golden eagle","bAccGen":"Northern goshawk","bColStr":"Speckled mousebird",
                        "bPogPus":"Red-fronted tinkerbird","bDryPub":"Downy woodpecker","bRisTri":"Black-legged kittiwake","bChrRid":"Black-headed gull",
                        "bGruAme":"Whooping crane","bGavSte":"Red-throated loon","bCalAnn":"Anna's hummingbird","bApuApu":"Common swift",
                        "bCucCan":"Common cuckoo","bCygOlo":"Mute swan","bAytFul":"Tufted duck","bAnaPla":"Pekin duck","bLagMut":"Rock ptarmigan",
                        "bGalGal":"White leghorn layer chicken","bRhePen":"Lesser rhea","rAllMis":"American alligator","rGopFla":"Mexican gopher tortoise",
                        "rMalTer":"Diamondback terrapin","rDerCor":"Leatherback sea turtle","rCheMyd":"Green sea turtle","rCarCar":"Loggerhead sea turtle",
                        "rEulEur":"European leaf-toed gecko","rThaEle":"Western terrestrial garter snake","rCanAsp":"Viper boa","rElgMul":"Southern alligator lizard",
                        "rRhiFlo":"Florida worm lizard","rZooViv":"European common lizard","rLacAgi":"Sand lizard","rPodRaf":"Aeolian wall lizard",
                        "aMicUni":"Tiny cayenne caecilian","aGeoSer":"Gaboon caecilian","aRhiBiv":"Two-lined caecilian","aBomBom":"European fire-bellied toad",
                        "aPelFus":"Common spadefoot toad","aRanTem":"European common frog","aBufBuf":"European toad","aHylSar":"Sardinian treefrog",
                        "aPseCor":"Southern corroboree frog"
                       }
            self.num_species=134
            self.numref_species=2
            self.ref_list =["pHomSap38","pHomSapT2T"]

        elif self.db=="t2t_16diploid":
            self.OrderList = ["hg38","hs1","GCA_018852605", "GCA_018852615",
                         "GCF_028858775","GCA_028858805",
                         "GCA_028858825","GCA_028858845",
                         "GCA_028885475", "GCA_028885495",
                         "GCA_028885625","GCA_028885525",
                         "GCA_028878055","GCA_028878085",
                         "GCA_028885655","GCA_028885685"]
            self.NameDict={"hg38":"hg38","hs1":"hs1",
                       "GCA_018852605":"hg002_pat",
                       "GCA_018852615":"hg002_mat",
                       "GCF_028858775": "Chimpanzee_hap1", 
                       "GCA_028858805":"Chimpanzee_hap2", 
                       "GCA_028858825":"P_Chimpanzee_pat",
                       "GCA_028858845":"P_Chimpanzee_mat",
                       "GCA_028885475":"Gorilla_pat",
                       "GCA_028885495":"Gorilla_mat",
                       "GCA_028885625":"Bornean_Orangutan_hap1",
                       "GCA_028885525":"Bornean_Orangutan_hap2",
                       "GCA_028878055":"S_Syndactylus_hap1",
                       "GCA_028878085":"S_Syndactylus_hap2",
                       "GCA_028885655":"Sumatran_Orangutan_hap1",
                       "GCA_028885685":"Sumatran_Orangutan_hap2"
                       }
            self.num_species=16
            self.numref_species=4
            self.ref_list =["hg38","hs1","GCA_018852605", "GCA_018852615"]
                           
        elif self.db=="t2t_hg38":
            self.OrderList = ["hg38","hs1","GCA_028858775","GCA_029289425","GCA_029281585","GCA_028885625","GCA_028885655","GCA_028878055"]
            self.NameDict={
                        "hg38":"hg38","hs1":"hs1","GCA_028858775": "chimpanzee", "GCA_029281585":"gorilla","GCA_029289425":"P_Chimpanzee",
                        "GCA_028885625":"Bornean_Orangutan","GCA_028878055":"S_syndactylus","GCA_028885655":"Sumatran_Orangutan"
                       }
            self.num_species=8
            self.numref_species=2
            self.ref_list =["hg38","hs1"]
        elif self.db=="t2tchimp":
            self.OrderList = ["GCA_028858775","hg38","hs1","GCA_029289425","GCA_029281585","GCA_028885625","GCA_028885655","GCA_028878055"]
            self.NameDict={
                        "hg38":"hg38","hs1":"hs1","GCA_028858775": "chimpanzee", "GCA_029281585":"gorilla","GCA_029289425":"P_Chimpanzee",
                        "GCA_028885625":"Bornean_Orangutan","GCA_028878055":"S_syndactylus","GCA_028885655":"Sumatran_Orangutan"
                       }
            self.num_species=8
            self.numref_species=1
            self.ref_list =["GCA_028858775"]
        elif self.db=="t2t_hs1":
            self.OrderList = ["hs1","hg38","GCA_028858775","GCA_029289425","GCA_029281585","GCA_028885625","GCA_028885655","GCA_028878055"]
            self.NameDict={
                        "hs1":"hs1","hg38":"hg38","GCA_028858775": "chimpanzee", "GCA_029281585":"gorilla","GCA_029289425":"P_Chimpanzee",
                        "GCA_028885625":"Bornean_Orangutan","GCA_028878055":"S_syndactylus","GCA_028885655":"Sumatran_Orangutan"
                       }
            self.num_species=8
            self.numref_species=2
            self.ref_list =["hs1","hg38"]
        elif self.db=="t2tgorilla":
            self.OrderList = ["GCA_029281585","hg38","hs1","GCA_028858775","GCA_029289425","GCA_028885625","GCA_028885655","GCA_028878055"]
            self.NameDict={
                        "hg38":"hg38","hs1":"hs1","GCA_028858775": "chimpanzee", "GCA_029281585":"gorilla","GCA_029289425":"P_Chimpanzee",
                        "GCA_028885625":"Bornean_Orangutan","GCA_028878055":"S_syndactylus","GCA_028885655":"Sumatran_Orangutan"
                       }
            self.num_species=8
            self.numref_species=1
            self.ref_list =["GCA_029281585"]
        elif self.db=="t2t_sorangutan":
            self.OrderList = ["GCA_028885655","hg38","hs1","GCA_028858775","GCA_029289425","GCA_029281585","GCA_028885625","GCA_028878055"]
            self.NameDict={
                        "hg38":"hg38","hs1":"hs1","GCA_028858775": "chimpanzee", "GCA_029281585":"gorilla","GCA_029289425":"P_Chimpanzee",
                        "GCA_028885625":"Bornean_Orangutan","GCA_028878055":"S_syndactylus","GCA_028885655":"Sumatran_Orangutan"
                       }
            self.num_species=8
            self.numref_species=1
            self.ref_list =["GCA_028885655"]
        elif self.db=="t2t_borangutan":
            self.OrderList = ["GCA_028885625","hg38","hs1","GCA_028858775","GCA_029289425","GCA_029281585","GCA_028885655","GCA_028878055"]
            self.NameDict={
                        "hg38":"hg38","hs1":"hs1","GCA_028858775": "chimpanzee", "GCA_029281585":"gorilla","GCA_029289425":"P_Chimpanzee",
                        "GCA_028885625":"Bornean_Orangutan","GCA_028878055":"S_syndactylus","GCA_028885655":"Sumatran_Orangutan"
                       }
            self.num_species=8
            self.numref_species=1
            self.ref_list =["GCA_028885625"]
        elif self.db=="t2t_siamang":
            self.OrderList = ["GCA_028878055","hg38","hs1","GCA_028858775","GCA_029289425","GCA_029281585","GCA_028885625","GCA_028885655"]
            self.NameDict={
                        "hg38":"hg38","hs1":"hs1","GCA_028858775": "chimpanzee", "GCA_029281585":"gorilla","GCA_029289425":"P_Chimpanzee",
                        "GCA_028885625":"Bornean_Orangutan","GCA_028878055":"S_syndactylus","GCA_028885655":"Sumatran_Orangutan"
                       }
            self.num_species=8
            self.numref_species=1
            self.ref_list =["GCA_028878055"]
        elif self.db=="gene_wide_ucsc":
            self.OrderList = ["hg38","hs1","Bonobo","Chimpanzee","Gorilla","Bornean_orangutan","Sumatran_orangutan","Siamang_gibbon","Common_marmoset","Sunda_slow_loris",
                         "Ring-tailed_lemur","Philippine_flying_lemur","American_pika","Grey_squirrel","Lesser_egyptian_jerboa","Deer_mouse","Southern_grasshopper_mouse",
                         "European_snow_vole","European_water_vole","Golden_spiny_mouse","Norway_rat","Nile_rat","European_woodmouse","House_mouse","Saddleback_dolphin",
                         "Bottlenose_dolphin","Long-finned_pilot_whale","White-beaked_dolphin","Killer_whale","Vaquita","Blainville's_beaked_whale","Pygmy_sperm_whale",
                         "North_atlantic_right_whale","Blue_whale","Common_minke_whale","Common_hippopotamus","Red_deer","Pig","Dromedary_camel","Arabian_horse",
                         "Southern_black_rhinoceros","Chinese_pangolin","Canada_lynx","Clouded_leopard","German_shepherd_dog","California_sea_lion","European_badger",
                         "Eurasian_otter","Stout","Greater_horseshoe_bat","Egyptian_fruit_bat","Pale_spear-nosed_bat","Common_vampire_bat","Pallas's_mastiff_bat",
                         "Kuhl's_pipistrelle","Daubenton's_bat","Greater_mouse-eared_bat","European_hedgehog","Etruscan_shrew","European_shrew","African_elephant",
                         "Asian_elephant","Nine-banded_armadillo","Linnaeus's_two-toed_sloth","Gray_short-tailed_opossum","Monito_del_montes","Tasmanian_devil",
                         "Common_brushtail_possum","Koala","Short-beaked_echidna","Platypus","Zebra_finch","Village_indigobird","Hawaiian_crow","New_caledonian_crow",
                         "Swamp_sparrow","Song_sparrow","Saltmarsh_sparrow","Brown-headed_cowbird","House_finch","Canary","White-throated_dipper","Swainson's_thrush",
                         "Black-capped_chickadee","Barn_swallow","Lanced-tailed_manakin","Budgerigar","Lesser_kestrel","Peregrine_falcon","Lanner_falcon",
                         "Gyrfalcon","Saker_falcon","European_golden_eagle","Northern_goshawk","Speckled_mousebird","Red-fronted_tinkerbird","Downy_woodpecker",
                         "Black-legged_kittiwake","Black-headed_gull","Whooping_crane","Red-throated_loon","Anna's_hummingbird","Common_swift","Common_cuckoo",
                         "Mute_swan","Tufted_duck","Pekin_duck","Rock_ptarmigan","White_leghorn_layer_chicken","Lesser_rhea","American_alligator","Mexican_gopher_tortoise",
                         "Diamondback_terrapin","Leatherback_sea_turtle","Green_sea_turtle","Loggerhead_sea_turtle","European_leaf-toed_gecko",
                         "Western_terrestrial_garter_snake","Viper_boa","Southern_alligator_lizard","Florida_worm_lizard","European_common_lizard",
                         "Sand_lizard","Aeolian_wall_lizard","Tiny_cayenne_caecilian","Gaboon_caecilian","Two-lined_caecilian","European_fire-bellied_toad",
                         "Common_spadefoot_toad","European_common_frog","European_toad","Sardinian_treefrog","Southern_corroboree_frog"]
            self.NameDict={
                        "pHomSap38":"hg38","pHomSapT2T":"hs1","pPanPan":"Bonobo","pPanTro":"Chimpanzee","pGorGor":"Gorilla","pPonPyg":"Bornean orangutan",
                        "pPonAbe":"Sumatran orangutan","pSymSyn":"Siamang gibbon","pCalJac":"Common marmoset","pNycCou":"Sunda slow loris","pLemCat":"Ring-tailed lemur",
                        "mCynVol":"Philippine flying lemur","mOchPri":"American pika","mSciCar":"Grey squirrel","mJacJac":"Lesser egyptian jerboa","mPerMan":"Deer mouse",
                        "mOnyTor":"Southern grasshopper mouse","mChiNiv":"European snow vole","mArvAmp":"European water vole","mAcoRus":"Golden spiny mouse",
                        "mRatNor":"Norway rat","mArvNil":"Nile rat","mApoSyl":"European woodmouse","mMusMus":"House mouse","mDelDel":"Saddleback dolphin",
                        "mTurTru":"Bottlenose dolphin","mGloMel":"Long-finned pilot whale","mLagAlb":"White-beaked dolphin","mOrcOrc":"Killer whale",
                        "mPhoSin":"Vaquita","mMesDen":"Blainville's beaked whale","mKogBre":"Pygmy sperm whale","mEubGla":"North atlantic right whale",
                        "mBalMus":"Blue whale","mBalAcu":"Common minke whale","mHipAmp":"Common hippopotamus","mBosTau":"Brahman cow","mCerEla":"Red deer",
                        "mSusScr":"Pig","mCamDro":"Dromedary camel","mEquCab":"Arabian horse","mDicBic":"Southern black rhinoceros","mManPen":"Chinese pangolin",
                        "mLynCan":"Canada lynx","mNeoNeb":"Clouded leopard","mCanLup":"German shepherd dog","mZalCal":"California sea lion","mMelMel":"European badger",
                        "mLutLut":"Eurasian otter","mMusErm":"Stout","mRhiFer":"Greater horseshoe bat","mRouAeg":"Egyptian fruit bat","mPhyDis":"Pale spear-nosed bat",
                        "mDesRot":"Common vampire bat","mMolMol":"Pallas's mastiff bat","mPipKuh":"Kuhl's pipistrelle","mMyoDau":"Daubenton's bat",
                        "mMyoMyo":"Greater mouse-eared bat","mEriEur":"European hedgehog","mSunEtr":"Etruscan shrew","mSorAra":"European shrew",
                        "mLoxAfr":"African elephant","mEleMax":"Asian elephant","mDasNov":"Nine-banded armadillo","mChoDid":"Linnaeus's two-toed sloth",
                        "mMonDom":"Gray short-tailed opossum","mDroGli":"Monito del montes","mSarHar":"Tasmanian devil","mTriVul":"Common brushtail possum",
                        "mPhaCin":"Koala","mTacAcu":"Short-beaked echidna","mOrnAna":"Platypus","bTaeGut":"Zebra finch","bVidCha":"Village indigobird",
                        "bCorHaw":"Hawaiian crow","bCorMon":"New caledonian crow","bMelGeo":"Swamp sparrow","bMelMel":"Song sparrow","bAmmCau":"Saltmarsh sparrow",
                        "bMolAte":"Brown-headed cowbird","bHaeMex":"House finch","bSerCan":"Canary","bCinCin":"White-throated dipper",
                        "bCatUst":"Swainson's thrush","bPoeAtr":"Black-capped chickadee","bHirRus":"Barn swallow","bChiLan":"Lanced-tailed manakin",
                        "bMelUnd":"Budgerigar","bFalNau":"Lesser kestrel","bFalPer":"Peregrine falcon","bFalBia":"Lanner falcon","bFalRus":"Gyrfalcon",
                        "bFalChe":"Saker falcon","bAquChr":"European golden eagle","bAccGen":"Northern goshawk","bColStr":"Speckled mousebird",
                        "bPogPus":"Red-fronted tinkerbird","bDryPub":"Downy woodpecker","bRisTri":"Black-legged kittiwake","bChrRid":"Black-headed gull",
                        "bGruAme":"Whooping crane","bGavSte":"Red-throated loon","bCalAnn":"Anna's hummingbird","bApuApu":"Common swift",
                        "bCucCan":"Common cuckoo","bCygOlo":"Mute swan","bAytFul":"Tufted duck","bAnaPla":"Pekin duck","bLagMut":"Rock ptarmigan",
                        "bGalGal":"White leghorn layer chicken","bRhePen":"Lesser rhea","rAllMis":"American alligator","rGopFla":"Mexican gopher tortoise",
                        "rMalTer":"Diamondback terrapin","rDerCor":"Leatherback sea turtle","rCheMyd":"Green sea turtle","rCarCar":"Loggerhead sea turtle",
                        "rEulEur":"European leaf-toed gecko","rThaEle":"Western terrestrial garter snake","rCanAsp":"Viper boa","rElgMul":"Southern alligator lizard",
                        "rRhiFlo":"Florida worm lizard","rZooViv":"European common lizard","rLacAgi":"Sand lizard","rPodRaf":"Aeolian wall lizard",
                        "aMicUni":"Tiny cayenne caecilian","aGeoSer":"Gaboon caecilian","aRhiBiv":"Two-lined caecilian","aBomBom":"European fire-bellied toad",
                        "aPelFus":"Common spadefoot toad","aRanTem":"European common frog","aBufBuf":"European toad","aHylSar":"Sardinian treefrog",
                        "aPseCor":"Southern corroboree frog"
                       }
            self.num_species=134
            self.numref_species=2
            self.ref_list =["hg38","hs1"]
        elif self.db=="gene_wide_ucsc_primate":
            self.OrderList = ["hg38","hs1","Bonobo","Chimpanzee","Gorilla","Bornean_orangutan","Sumatran_orangutan","Siamang_gibbon","Common_marmoset","Sunda_slow_loris",
                         "Ring-tailed_lemur","Philippine_flying_lemur"]
            self.NameDict={
                        "pHomSap38":"hg38","pHomSapT2T":"hs1","pPanPan":"Bonobo","pPanTro":"Chimpanzee","pGorGor":"Gorilla","pPonPyg":"Bornean orangutan",
                        "pPonAbe":"Sumatran orangutan","pSymSyn":"Siamang gibbon","pCalJac":"Common marmoset","pNycCou":"Sunda slow loris","pLemCat":"Ring-tailed lemur",
                        "mCynVol":"Philippine flying lemur","mOchPri":"American pika","mSciCar":"Grey squirrel","mJacJac":"Lesser egyptian jerboa","mPerMan":"Deer mouse",
                        "mOnyTor":"Southern grasshopper mouse","mChiNiv":"European snow vole","mArvAmp":"European water vole","mAcoRus":"Golden spiny mouse",
                        "mRatNor":"Norway rat","mArvNil":"Nile rat","mApoSyl":"European woodmouse","mMusMus":"House mouse","mDelDel":"Saddleback dolphin",
                        "mTurTru":"Bottlenose dolphin","mGloMel":"Long-finned pilot whale","mLagAlb":"White-beaked dolphin","mOrcOrc":"Killer whale",
                        "mPhoSin":"Vaquita","mMesDen":"Blainville's beaked whale","mKogBre":"Pygmy sperm whale","mEubGla":"North atlantic right whale",
                        "mBalMus":"Blue whale","mBalAcu":"Common minke whale","mHipAmp":"Common hippopotamus","mBosTau":"Brahman cow","mCerEla":"Red deer",
                        "mSusScr":"Pig","mCamDro":"Dromedary camel","mEquCab":"Arabian horse","mDicBic":"Southern black rhinoceros","mManPen":"Chinese pangolin",
                        "mLynCan":"Canada lynx","mNeoNeb":"Clouded leopard","mCanLup":"German shepherd dog","mZalCal":"California sea lion","mMelMel":"European badger",
                        "mLutLut":"Eurasian otter","mMusErm":"Stout","mRhiFer":"Greater horseshoe bat","mRouAeg":"Egyptian fruit bat","mPhyDis":"Pale spear-nosed bat",
                        "mDesRot":"Common vampire bat","mMolMol":"Pallas's mastiff bat","mPipKuh":"Kuhl's pipistrelle","mMyoDau":"Daubenton's bat",
                        "mMyoMyo":"Greater mouse-eared bat","mEriEur":"European hedgehog","mSunEtr":"Etruscan shrew","mSorAra":"European shrew",
                        "mLoxAfr":"African elephant","mEleMax":"Asian elephant","mDasNov":"Nine-banded armadillo","mChoDid":"Linnaeus's two-toed sloth",
                        "mMonDom":"Gray short-tailed opossum","mDroGli":"Monito del montes","mSarHar":"Tasmanian devil","mTriVul":"Common brushtail possum",
                        "mPhaCin":"Koala","mTacAcu":"Short-beaked echidna","mOrnAna":"Platypus","bTaeGut":"Zebra finch","bVidCha":"Village indigobird",
                        "bCorHaw":"Hawaiian crow","bCorMon":"New caledonian crow","bMelGeo":"Swamp sparrow","bMelMel":"Song sparrow","bAmmCau":"Saltmarsh sparrow",
                        "bMolAte":"Brown-headed cowbird","bHaeMex":"House finch","bSerCan":"Canary","bCinCin":"White-throated dipper",
                        "bCatUst":"Swainson's thrush","bPoeAtr":"Black-capped chickadee","bHirRus":"Barn swallow","bChiLan":"Lanced-tailed manakin",
                        "bMelUnd":"Budgerigar","bFalNau":"Lesser kestrel","bFalPer":"Peregrine falcon","bFalBia":"Lanner falcon","bFalRus":"Gyrfalcon",
                        "bFalChe":"Saker falcon","bAquChr":"European golden eagle","bAccGen":"Northern goshawk","bColStr":"Speckled mousebird",
                        "bPogPus":"Red-fronted tinkerbird","bDryPub":"Downy woodpecker","bRisTri":"Black-legged kittiwake","bChrRid":"Black-headed gull",
                        "bGruAme":"Whooping crane","bGavSte":"Red-throated loon","bCalAnn":"Anna's hummingbird","bApuApu":"Common swift",
                        "bCucCan":"Common cuckoo","bCygOlo":"Mute swan","bAytFul":"Tufted duck","bAnaPla":"Pekin duck","bLagMut":"Rock ptarmigan",
                        "bGalGal":"White leghorn layer chicken","bRhePen":"Lesser rhea","rAllMis":"American alligator","rGopFla":"Mexican gopher tortoise",
                        "rMalTer":"Diamondback terrapin","rDerCor":"Leatherback sea turtle","rCheMyd":"Green sea turtle","rCarCar":"Loggerhead sea turtle",
                        "rEulEur":"European leaf-toed gecko","rThaEle":"Western terrestrial garter snake","rCanAsp":"Viper boa","rElgMul":"Southern alligator lizard",
                        "rRhiFlo":"Florida worm lizard","rZooViv":"European common lizard","rLacAgi":"Sand lizard","rPodRaf":"Aeolian wall lizard",
                        "aMicUni":"Tiny cayenne caecilian","aGeoSer":"Gaboon caecilian","aRhiBiv":"Two-lined caecilian","aBomBom":"European fire-bellied toad",
                        "aPelFus":"Common spadefoot toad","aRanTem":"European common frog","aBufBuf":"European toad","aHylSar":"Sardinian treefrog",
                        "aPseCor":"Southern corroboree frog"
                       }
            self.num_species=12
            self.numref_species=2
            self.ref_list =["hg38","hs1"]
        elif self.db=="gene_wide_ucsc_bird":
            self.OrderList = ["Zebra_finch","Village_indigobird","Hawaiian_crow","New_caledonian_crow",
                         "Swamp_sparrow","Song_sparrow","Saltmarsh_sparrow","Brown-headed_cowbird","House_finch","Canary","White-throated_dipper","Swainson's_thrush",
                         "Black-capped_chickadee","Barn_swallow","Budgerigar","Anna's_hummingbird","Lanced-tailed_manakin","Lesser_kestrel","Peregrine_falcon","Lanner_falcon",
                         "Gyrfalcon","Saker_falcon","European_golden_eagle","Northern_goshawk","Speckled_mousebird","Red-fronted_tinkerbird","Downy_woodpecker",
                         "Black-legged_kittiwake","Black-headed_gull","Whooping_crane","Red-throated_loon","Common_swift","Common_cuckoo",
                         "Mute_swan","Tufted_duck","Pekin_duck","Rock_ptarmigan","White_leghorn_layer_chicken","Lesser_rhea"]
            self.NameDict={
                        "pHomSap38":"hg38","pHomSapT2T":"hs1","pPanPan":"Bonobo","pPanTro":"Chimpanzee","pGorGor":"Gorilla","pPonPyg":"Bornean orangutan",
                        "pPonAbe":"Sumatran orangutan","pSymSyn":"Siamang gibbon","pCalJac":"Common marmoset","pNycCou":"Sunda slow loris","pLemCat":"Ring-tailed lemur",
                        "mCynVol":"Philippine flying lemur","mOchPri":"American pika","mSciCar":"Grey squirrel","mJacJac":"Lesser egyptian jerboa","mPerMan":"Deer mouse",
                        "mOnyTor":"Southern grasshopper mouse","mChiNiv":"European snow vole","mArvAmp":"European water vole","mAcoRus":"Golden spiny mouse",
                        "mRatNor":"Norway rat","mArvNil":"Nile rat","mApoSyl":"European woodmouse","mMusMus":"House mouse","mDelDel":"Saddleback dolphin",
                        "mTurTru":"Bottlenose dolphin","mGloMel":"Long-finned pilot whale","mLagAlb":"White-beaked dolphin","mOrcOrc":"Killer whale",
                        "mPhoSin":"Vaquita","mMesDen":"Blainville's beaked whale","mKogBre":"Pygmy sperm whale","mEubGla":"North atlantic right whale",
                        "mBalMus":"Blue whale","mBalAcu":"Common minke whale","mHipAmp":"Common hippopotamus","mBosTau":"Brahman cow","mCerEla":"Red deer",
                        "mSusScr":"Pig","mCamDro":"Dromedary camel","mEquCab":"Arabian horse","mDicBic":"Southern black rhinoceros","mManPen":"Chinese pangolin",
                        "mLynCan":"Canada lynx","mNeoNeb":"Clouded leopard","mCanLup":"German shepherd dog","mZalCal":"California sea lion","mMelMel":"European badger",
                        "mLutLut":"Eurasian otter","mMusErm":"Stout","mRhiFer":"Greater horseshoe bat","mRouAeg":"Egyptian fruit bat","mPhyDis":"Pale spear-nosed bat",
                        "mDesRot":"Common vampire bat","mMolMol":"Pallas's mastiff bat","mPipKuh":"Kuhl's pipistrelle","mMyoDau":"Daubenton's bat",
                        "mMyoMyo":"Greater mouse-eared bat","mEriEur":"European hedgehog","mSunEtr":"Etruscan shrew","mSorAra":"European shrew",
                        "mLoxAfr":"African elephant","mEleMax":"Asian elephant","mDasNov":"Nine-banded armadillo","mChoDid":"Linnaeus's two-toed sloth",
                        "mMonDom":"Gray short-tailed opossum","mDroGli":"Monito del montes","mSarHar":"Tasmanian devil","mTriVul":"Common brushtail possum",
                        "mPhaCin":"Koala","mTacAcu":"Short-beaked echidna","mOrnAna":"Platypus","bTaeGut":"Zebra finch","bVidCha":"Village indigobird",
                        "bCorHaw":"Hawaiian crow","bCorMon":"New caledonian crow","bMelGeo":"Swamp sparrow","bMelMel":"Song sparrow","bAmmCau":"Saltmarsh sparrow",
                        "bMolAte":"Brown-headed cowbird","bHaeMex":"House finch","bSerCan":"Canary","bCinCin":"White-throated dipper",
                        "bCatUst":"Swainson's thrush","bPoeAtr":"Black-capped chickadee","bHirRus":"Barn swallow","bChiLan":"Lanced-tailed manakin",
                        "bMelUnd":"Budgerigar","bFalNau":"Lesser kestrel","bFalPer":"Peregrine falcon","bFalBia":"Lanner falcon","bFalRus":"Gyrfalcon",
                        "bFalChe":"Saker falcon","bAquChr":"European golden eagle","bAccGen":"Northern goshawk","bColStr":"Speckled mousebird",
                        "bPogPus":"Red-fronted tinkerbird","bDryPub":"Downy woodpecker","bRisTri":"Black-legged kittiwake","bChrRid":"Black-headed gull",
                        "bGruAme":"Whooping crane","bGavSte":"Red-throated loon","bCalAnn":"Anna's hummingbird","bApuApu":"Common swift",
                        "bCucCan":"Common cuckoo","bCygOlo":"Mute swan","bAytFul":"Tufted duck","bAnaPla":"Pekin duck","bLagMut":"Rock ptarmigan",
                        "bGalGal":"White leghorn layer chicken","bRhePen":"Lesser rhea","rAllMis":"American alligator","rGopFla":"Mexican gopher tortoise",
                        "rMalTer":"Diamondback terrapin","rDerCor":"Leatherback sea turtle","rCheMyd":"Green sea turtle","rCarCar":"Loggerhead sea turtle",
                        "rEulEur":"European leaf-toed gecko","rThaEle":"Western terrestrial garter snake","rCanAsp":"Viper boa","rElgMul":"Southern alligator lizard",
                        "rRhiFlo":"Florida worm lizard","rZooViv":"European common lizard","rLacAgi":"Sand lizard","rPodRaf":"Aeolian wall lizard",
                        "aMicUni":"Tiny cayenne caecilian","aGeoSer":"Gaboon caecilian","aRhiBiv":"Two-lined caecilian","aBomBom":"European fire-bellied toad",
                        "aPelFus":"Common spadefoot toad","aRanTem":"European common frog","aBufBuf":"European toad","aHylSar":"Sardinian treefrog",
                        "aPseCor":"Southern corroboree frog"
                       }
            self.num_species=39
            self.numref_species=16
            self.ref_list =["Zebra_finch","Village_indigobird","Hawaiian_crow","New_caledonian_crow",
                         "Swamp_sparrow","Song_sparrow","Saltmarsh_sparrow","Brown-headed_cowbird","House_finch","Canary","White-throated_dipper","Swainson's_thrush",
                         "Black-capped_chickadee","Barn_swallow","Budgerigar","Anna's_hummingbird"]
    def setNameList(self, convertlist):
        self.nameList=[]
        for index in range(0,len(convertlist)):
            #print(" appending "+self.NameDict.get(convertlist[index])+" at index "+ str(index)+"\n")
            self.nameList.append(self.NameDict.get(convertlist[index]))
        return self.nameList
    
    def getName(self, geneName):
        CommonName=self.NameDict[geneName]
        return CommonName
    def getRefList(self):
        RefList=self.ref_list
        return RefList
    def getNumRefSpecies(self):
        numRef=self.numref_species
        return numRef
    