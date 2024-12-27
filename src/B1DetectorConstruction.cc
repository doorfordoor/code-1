//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class
//*
#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

G4VPhysicalVolume* B1DetectorConstruction::Construct()
//Constuct不可改动，返回值G4VPhysicalVolume
{  
G4NistManager* nist = G4NistManager::Instance();
  G4double env_sizeX = 40*mm, env_sizeY = 40*mm, env_sizeZ = 10*mm;
  G4double world_sizeXY = 100*mm ,world_sizeZ = 50*mm;


  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* copper = nist->FindOrBuildMaterial("G4_Cu");   
  G4double p, a; 
  G4Element* ele_H  = new G4Element("Hydrogen","H" ,p=1. , a =1.00794*g/mole);
  G4Element* ele_C  = new G4Element("Carbon"  ,"C" ,p=6. , a =12.01 *g/mole);
  G4Element* ele_O  = new G4Element("Oxygen"  ,"O" ,p=8. , a =16.00 *g/mole);
  G4Element* ele_N  = new G4Element("Nitrogen","N" ,p=7. , a =14.01 *g/mole); 
  G4Element* ele_Si = new G4Element("Silicon" ,"Si",p=14., a =28.086*g/mole); 
  G4Material* Filled= new G4Material("Filled",1.20*g/cm3, 3);
  Filled->AddElement(ele_C, 11);
  Filled->AddElement(ele_H, 12);
  Filled->AddElement(ele_O, 3);
  //1.20*g/cm3, 3
  G4Material* PolyimideResin = new G4Material("PolyimideResin", 1.20*g/cm3, 4);
  PolyimideResin->AddElement(ele_C, 35);
  PolyimideResin->AddElement(ele_H, 28);
  PolyimideResin->AddElement(ele_N, 2);
  PolyimideResin->AddElement(ele_O, 7);
  //G4Material* SiliconDioxide = nist -> FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* SiliconDioxide = new G4Material("SiliconDioxide", 2.20*g/cm3, 2);
  SiliconDioxide->AddElement(ele_Si, 1);
  SiliconDioxide->AddElement(ele_O, 2);
  G4Material* EpoxyResin = new G4Material("EpoxyResin", 1.70*g/cm3, 2);
  EpoxyResin->AddMaterial(Filled, 50.0*perCent);
  EpoxyResin->AddMaterial(SiliconDioxide, 50.0*perCent);
  G4Material* Silicon = nist -> FindOrBuildMaterial("G4_Si");
  //定义材料

  //定义世界层
  G4Box* solidWorld = 
    new G4Box("World",
    0.5*world_sizeXY,0.5*world_sizeXY,0.5*world_sizeZ);
  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,
                        air,
                        "World");
  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,
                      G4ThreeVector(),
                      logicWorld,
                      "World",
                      0,
                      false,
                      0,
                      true);
 
  //定义环境层                     
  G4Box* solidEnv =
    new G4Box("Env",
    0.5*env_sizeX,0.5*env_sizeY,0.5*env_sizeZ);
  G4LogicalVolume* logicEnv =
    new G4LogicalVolume(solidEnv,
                        air,
                        "Env");
    new G4PVPlacement(0,
                      G4ThreeVector(),
                      logicEnv,
                      "Env",
                      logicWorld,
                      false,
                      0,
                      true);
  
  //UBM_B_5层,type:Metal,Material:Air
  G4double chip_sizeX = 10.77*mm, chip_sizeY = 9.90*mm;
  G4double UBM_B_5_sizeZ = 0.004*mm;
  G4Box* solidUBM_B_5 =
    new G4Box("UBM_B_5",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*UBM_B_5_sizeZ);
  G4LogicalVolume* logicUBM_B_5 =
    new G4LogicalVolume(solidUBM_B_5,
                        air,
                        "UBM_B_5");
  G4ThreeVector posUBM_B_5 = G4ThreeVector(0 ,0 ,0.5*UBM_B_5_sizeZ);
    new G4PVPlacement(0,
                      posUBM_B_5,
                      logicUBM_B_5,
                      "UBM_B_5",
                      logicEnv,
                      false,
                      0,
                      true);

// //UBM_B_5层的焊接TSVOUTERBGAPADa
// G4double TSVOUTERBGAPAD_Ra =0.135*mm;
// G4Tubs* solidTSVOUTERBGAPADa =
//     new G4Tubs("TSVOUTERBGAPADa",
//     0, 
//     TSVOUTERBGAPAD_Ra, 
//     0.5*UBM_B_5_sizeZ, 
//     0.*deg, 
//     360.*deg);
// G4LogicalVolume* logicTSVOUTERBGAPADa =
//     new G4LogicalVolume(solidTSVOUTERBGAPADa,
//                         copper,
//                         "TSVOUTERBGAPADa");
// G4int nTSVOUTERBGAPADa = 0;
// for (G4double x = -4.875*mm ; x <= 4.875*mm ; x = x + 0.65*mm){
//   for (G4double y = -4.55*mm ; y <= 4.55*mm ; y = y + 0.65*mm){
//     G4ThreeVector posTSVOUTERBGAPADa = G4ThreeVector(x, y, 0);
//     new G4PVPlacement(0,
//                       posTSVOUTERBGAPADa,
//                       logicTSVOUTERBGAPADa,
//                       "TSVOUTERBGAPADa",
//                       logicUBM_B_5,
//                       false,
//                       nTSVOUTERBGAPADa++,
//                       true);
//   }
// };

//BP2_5层,type:Dielectric,Material:FR_4
G4double BP2_5High = 0.2072*mm;
G4Box* solidBP2_5 =
    new G4Box("BP2_5",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BP2_5High-UBM_B_5_sizeZ));
G4LogicalVolume* logicBP2_5 =
    new G4LogicalVolume(solidBP2_5,
                        EpoxyResin,
                        "BP2_5");
G4ThreeVector posBP2_5 = G4ThreeVector(0,0,0.5*(BP2_5High+UBM_B_5_sizeZ));
    new G4PVPlacement(0,
                      posBP2_5,
                      logicBP2_5,
                      "BP2_5",
                      logicEnv,
                      false,
                      0,
                      true);

// //BP2_5层的焊接TSVOUTERBGAPADb
// G4double TSVOUTERBGAPAD_Rb =0.135*mm;
// G4Tubs* solidTSVOUTERBGAPADb =
//     new G4Tubs("TSVOUTERBGAPADb",
//     0, 
//     TSVOUTERBGAPAD_Rb, 
//     0.5*(BP2_5High-UBM_B_5_sizeZ), 
//     0.*deg, 
//     360.*deg);
// G4LogicalVolume* logicTSVOUTERBGAPADb =
//     new G4LogicalVolume(solidTSVOUTERBGAPADb,
//                         copper,
//                         "TSVOUTERBGAPADb");
// G4int nTSVOUTERBGAPADb = 0;
// for (G4double x = -4.875*mm ; x <= 4.875*mm ; x = x + 0.65*mm){
//   for (G4double y = -4.55*mm ; y <= 4.55*mm ; y = y + 0.65*mm){
//     G4ThreeVector posTSVOUTERBGAPADb = G4ThreeVector(x, y, 0);
//     new G4PVPlacement(0,
//                       posTSVOUTERBGAPADb,
//                       logicTSVOUTERBGAPADb,
//                       "TSVOUTERBGAPADb",
//                       logicBP2_5,
//                       false,
//                       nTSVOUTERBGAPADb++,
//                       true);
//   }
// };

//BM2_5层,type:Metal,Material:Filled
G4double BM2_5High = 0.23768*mm;
//G4double BP2_5High = 0.2072*mm;
G4Box* solidBM2_5 =
    new G4Box("BM2_5",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BM2_5High-BP2_5High));

G4LogicalVolume* logicBM2_5 =
    new G4LogicalVolume(solidBM2_5,
                        Filled,
                        "BM2_5");
G4ThreeVector posBM2_5 = G4ThreeVector(0,0,0.5*(BM2_5High+BP2_5High));
    new G4PVPlacement(0,
                      posBM2_5,
                      logicBM2_5,
                      "BM2_5",
                      logicEnv,
                      false,
                      0,
                      true);

// //BM2_5层的焊接TSVOUTERBGAPADc
// G4double TSVOUTERBGAPAD_Rc =0.135*mm;
// G4Tubs* solidTSVOUTERBGAPADc =
//     new G4Tubs("TSVOUTERBGAPADc",
//     0, 
//     TSVOUTERBGAPAD_Rc, 
//     0.5*(BM2_5High-BP2_5High), 
//     0.*deg, 
//     360.*deg);
// G4LogicalVolume* logicTSVOUTERBGAPADc =
//     new G4LogicalVolume(solidTSVOUTERBGAPADc,
//                         copper,
//                         "TSVOUTERBGAPADc");
// G4int nTSVOUTERBGAPADc = 0;
// for (G4double x = -4.875*mm ; x <= 4.875*mm ; x = x + 0.65*mm){
//   for (G4double y = -4.55*mm ; y <= 4.55*mm ; y = y + 0.65*mm){
//     G4ThreeVector posTSVOUTERBGAPADc = G4ThreeVector(x, y, 0);
//     new G4PVPlacement(0,
//                       posTSVOUTERBGAPADc,
//                       logicTSVOUTERBGAPADc,
//                       "TSVOUTERBGAPADc",
//                       logicBM2_5,
//                       false,
//                       nTSVOUTERBGAPADc++,
//                       true);
//   }
// };

//TSVOUTERBGAPAD层
//BM2_5High = 0.23768*mm;
G4double TSVOUTERBGAPAD_R =0.135*mm;
G4Tubs* solidTSVOUTERBGAPAD =
    new G4Tubs("TSVOUTERBGAPAD",
    0, 
    TSVOUTERBGAPAD_R, 
    0.5*BM2_5High, 
    0.*deg, 
    360.*deg);
G4LogicalVolume* logicTSVOUTERBGAPAD =
    new G4LogicalVolume(solidTSVOUTERBGAPAD,
                        copper,
                        "TSVOUTERBGAPAD");
G4int nTSVOUTERBGAPAD = 0;
for (G4double x = -4.875*mm ; x <= 4.875*mm ; x = x + 0.65*mm){
  for (G4double y = -4.55*mm ; y <= 4.55*mm ; y = y + 0.65*mm){
    G4ThreeVector posTSVOUTERBGAPAD = G4ThreeVector(x, y, 0.5*BM2_5High);
    new G4PVPlacement(0,
                      posTSVOUTERBGAPAD,
                      logicTSVOUTERBGAPAD,
                      "TSVOUTERBGAPAD",
                      logicEnv,
                      false,
                      nTSVOUTERBGAPAD++,
                      true);
  }
};

//BP1_5层,type:Dielectric,Material:FR_4
G4double BP1_5High = 0.24768*mm;
G4Box* solidBP1_5 =
    new G4Box("BP1_5",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BP1_5High-BM2_5High));
G4LogicalVolume* logicBP1_5 =
  new G4LogicalVolume(solidBP1_5,
                      EpoxyResin,
                     "BP1_5");
G4ThreeVector posBP1_5 = G4ThreeVector(0,0,0.5*(BP1_5High+BM2_5High));
    new G4PVPlacement(0,
                      posBP1_5,
                      logicBP1_5,
                      "BP1_5",
                      logicEnv,
                      false,
                      0,
                      true);

//BM1_5层,type:Metal,Material:Filled
G4double BM1_5High = 0.24968*mm;
G4Box* solidBM1_5 =
    new G4Box("BM1_5",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BM1_5High-BP1_5High));
    G4LogicalVolume* logicBM1_5 =
    new G4LogicalVolume(solidBM1_5,
                        Filled,
                        "BM1_5");
G4ThreeVector posBM1_5 = G4ThreeVector(0,0,0.5*(BM1_5High+BP1_5High));
    new G4PVPlacement(0,
                      posBM1_5,
                      logicBM1_5,
                      "BM1_5",
                      logicEnv,
                      false,
                      0,
                      true);

//TSV-BBVIA30-BM1_5_BM2_5,copper,r=0.015mm 
//BM1_5High = 0.24968*mm;
G4double TSV_BBVIA30_BM1_5_BM2_5R =0.015*mm;
G4Tubs* solidTSV_BBVIA30_BM1_5_BM2_5 =
    new G4Tubs("TSV_BBVIA30_BM1_5_BM2_5",
    0,
    TSV_BBVIA30_BM1_5_BM2_5R,
    0.5*(BM1_5High-BM2_5High),
    0.*deg,
    360.*deg); 
G4LogicalVolume* logicTSV_BBVIA30_BM1_5_BM2_5 =
    new G4LogicalVolume(solidTSV_BBVIA30_BM1_5_BM2_5,
                        copper,
                        "TSV_BBVIA30_BM1_5_BM2_5");
//放置527个TSV_BBVIA30，且不能和TSVOUTERBGAPAD以及其他TSV_BBVIA30重叠
// for (int i = 0; i < 527; ++i) {
//     G4double x = (G4UniformRand() - 0.5) * (10.57*mm - 2 * TSV_BBVIA30_BM1_5_BM2_5R);
//     G4double y = (G4UniformRand() - 0.5) * (9.7*mm - 2 * TSV_BBVIA30_BM1_5_BM2_5R);
//     G4double z = 0.5 * (BM1_5High + BM2_5High);
//     G4ThreeVector posTSV_BBVIA30_BM1_5_BM2_5(x, y, z);
//     while (OverlapWithOtherTSV_BBVIA30_BM1_5_BM2_5(posTSV_BBVIA30_BM1_5_BM2_5, TSV_BBVIA30_BM1_5_BM2_5R)) {
//      G4double x = (G4UniformRand() - 0.5) * (10.57*mm - 2 * TSV_BBVIA30_BM1_5_BM2_5R);
//      G4double y = (G4UniformRand() - 0.5) * (9.7*mm - 2 * TSV_BBVIA30_BM1_5_BM2_5R);
//         posTSV_BBVIA30_BM1_5_BM2_5 = G4ThreeVector(x, y, z);
//     }
//     new G4PVPlacement(0,
//                       posTSV_BBVIA30_BM1_5_BM2_5,
//                       logicTSV_BBVIA30_BM1_5_BM2_5,
//                       "TSV_BBVIA30_BM1_5_BM2_5",
//                       logicEnv,
//                       false,
//                       i,
//                       true);
// }
// 计算每行可以放置的TSV数量

// 计算TSV之间的水平和垂直间距
G4double spacingX = (10.57 * mm - 2 * TSV_BBVIA30_BM1_5_BM2_5R) / 30;
G4double spacingY = (9.7 * mm - 2 * TSV_BBVIA30_BM1_5_BM2_5R) / 16;

for (int i = 0; i < 526; ++i) {
    // 计算当前TSV在网格中的行和列
    int row = i / 31;
    int col = i % 31;

    // 计算当前TSV的位置
    G4double x = -0.5 * (10.57 * mm - 2 * TSV_BBVIA30_BM1_5_BM2_5R) + col * spacingX;
    G4double y = -0.5 * (9.7 * mm - 2 * TSV_BBVIA30_BM1_5_BM2_5R) + row * spacingY;
    G4double z = 0.5 * (BM1_5High + BM2_5High);
    G4ThreeVector posTSV_BBVIA30_BM1_5_BM2_5(x, y, z);

    // 检查是否与其他TSV重叠
    // if (OverlapWithOtherTSV_BBVIA30_BM1_5_BM2_5(posTSV_BBVIA30_BM1_5_BM2_5, TSV_BBVIA30_BM1_5_BM2_5R)) {
    //     // 如果重叠，打印错误信息并停止程序
    //     G4cerr << "TSV at position (" << x << ", " << y << ", " << z << ") overlaps with another TSV!" << G4endl;
    //     exit(1);
    // }

    // 创建并放置TSV
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_5_BM2_5,
                      logicTSV_BBVIA30_BM1_5_BM2_5,
                      "TSV_BBVIA30_BM1_5_BM2_5",
                      logicEnv,
                      false,
                      i,
                      true);
}




//BPI_CORE_5层,type:Dielectric,Material:FR_4                      
G4double BPI_CORE_5High = 0.25768*mm;
G4Box* solidBPI_CORE_5 =
    new G4Box("BPI_CORE_5",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BPI_CORE_5High-BM1_5High));
G4LogicalVolume* logicBPI_CORE_5 =
    new G4LogicalVolume(solidBPI_CORE_5,
                        EpoxyResin,
                        "BPI_CORE_5");
G4ThreeVector posBPI_CORE_5 = G4ThreeVector(0,0,0.5*(BPI_CORE_5High+BM1_5High));
    new G4PVPlacement(0,
                      posBPI_CORE_5,
                      logicBPI_CORE_5,
                      "BPI_CORE_5",
                      logicEnv,
                      false,
                      0,
                      true);

//CORE2_5层，type:Dielectric,Material:FR_4
G4double CORE2_5High = 0.28768*mm;
G4Box* solidCORE2_5 =
    new G4Box("CORE2_5",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(CORE2_5High-BPI_CORE_5High));
G4LogicalVolume* logicCORE2_5 =
    new G4LogicalVolume(solidCORE2_5,
                        EpoxyResin,
                        "CORE2_5");
G4ThreeVector posCORE2_5 = G4ThreeVector(0,0,0.5*(CORE2_5High+BPI_CORE_5High));
    new G4PVPlacement(0,
                      posCORE2_5,
                      logicCORE2_5,
                      "CORE2_5",
                      logicEnv,
                      false,
                      0,
                      true);

//CORE1_5层，type:Dielectric,Material:FR_4
G4double CORE1_5High = 0.40768*mm;
G4Box* solidCORE1_5 =
    new G4Box("CORE1_5",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(CORE1_5High-CORE2_5High));
G4LogicalVolume* logicCORE1_5 =
    new G4LogicalVolume(solidCORE1_5,
                        EpoxyResin,
                        "CORE1_5");
G4ThreeVector posCORE1_5 = G4ThreeVector(0,0,0.5*(CORE1_5High+CORE2_5High));
    new G4PVPlacement(0,
                      posCORE1_5,
                      logicCORE1_5,
                      "CORE1_5",
                      logicEnv,
                      false,
                      0,
                      true);                            

//FPI_CORE_5层，type:Dielectric,Material:FR_4
G4double FPI_CORE_5High = 0.41568*mm;
G4Box* solidFPI_CORE_5 =
    new G4Box("FPI_CORE_5",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FPI_CORE_5High-CORE1_5High));
G4LogicalVolume* logicFPI_CORE_5 =
    new G4LogicalVolume(solidFPI_CORE_5,
                        EpoxyResin,
                        "FPI_CORE_5");
G4ThreeVector posFPI_CORE_5 = G4ThreeVector(0,0,0.5*(FPI_CORE_5High+CORE1_5High));
    new G4PVPlacement(0,
                      posFPI_CORE_5,
                      logicFPI_CORE_5,
                      "FPI_CORE_5",
                      logicEnv,
                      false,
                      0,
                      true);

//FM1_5层，type:Metal,Material:Filled
G4double FM1_5High = 0.41968*mm;
G4Box* solidFM1_5 =
    new G4Box("FM1_5",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FM1_5High-FPI_CORE_5High));
G4LogicalVolume* logicFM1_5 =
    new G4LogicalVolume(solidFM1_5,
                        Filled,
                        "FM1_5");
G4ThreeVector posFM1_5 = G4ThreeVector(0,0,0.5*(FM1_5High+FPI_CORE_5High));
    new G4PVPlacement(0,
                      posFM1_5,
                      logicFM1_5,
                      "FM1_5",
                      logicEnv,
                      false,
                      0,
                      true);

//焊盘：TSV-BBVIA-30-FM1_5_BM1_5,copper,679个，r=0.01mm
G4double TSV_BBVIA30_FM1_5_BM1_5R = 0.01*mm;
G4Tubs* solidTSV_BBVIA30_FM1_5_BM1_5 =
    new G4Tubs("TSV_BBVIA30_FM1_5_BM1_5",
    0,
    TSV_BBVIA30_FM1_5_BM1_5R,
    0.5*(FM1_5High-BM1_5High),
    0.*deg,
    360.*deg);
G4LogicalVolume* logicTSV_BBVIA30_FM1_5_BM1_5 =
    new G4LogicalVolume(solidTSV_BBVIA30_FM1_5_BM1_5,
                        copper,
                        "TSV_BBVIA30_FM1_5_BM1_5");
G4int numTSVs1 = 192;
G4double boxXMin1 = -5.285 * mm + TSV_BBVIA30_FM1_5_BM1_5R;
G4double boxXMax1 = 5.285 * mm - TSV_BBVIA30_FM1_5_BM1_5R;
G4double boxYMin1 = 4.065 * mm + 5*TSV_BBVIA30_FM1_5_BM1_5R;
G4double boxYMax1 = 4.85 * mm - TSV_BBVIA30_FM1_5_BM1_5R;  
G4double xStep1 = (boxXMax1 - boxXMin1) / 47;// 48个TSV_BBVIA30_FM1_5_BM1_5
G4double yStep1 = (boxYMax1 - boxYMin1) / 3; //4列TSV_BBVIA30_FM1_5_BM1_5
for (int i = 0; i < numTSVs1; ++i) {
    G4double x = boxXMin1 + (i % 48) * xStep1;
    G4double y = boxYMin1 + (i / 48) * yStep1;
    G4double z = 0.5 * (FM1_5High + BM1_5High);
    G4ThreeVector posTSV_BBVIA30_FM1_5_BM1_5(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA30_FM1_5_BM1_5,
                      logicTSV_BBVIA30_FM1_5_BM1_5,
                      "TSV_BBVIA30_FM1_5_BM1_5",
                      logicEnv,
                      false,
                      i,
                      true);
}
G4int numTSVs2 = 192;    
G4double boxXMin2 = -5.285 * mm + TSV_BBVIA30_FM1_5_BM1_5R;
G4double boxXMax2 = 5.285 * mm - TSV_BBVIA30_FM1_5_BM1_5R;
G4double boxYMin2 = -4.85 * mm + TSV_BBVIA30_FM1_5_BM1_5R;  
G4double boxYMax2 = -4.065 * mm - 5*TSV_BBVIA30_FM1_5_BM1_5R;               
G4double xStep2 = (boxXMax2 - boxXMin2) / 47;
G4double yStep2 = (boxYMax2 - boxYMin2) / 3;
for (int i = 0; i < numTSVs2; ++i) {
    G4double x = boxXMin2 + (i % 48) * xStep2;
    G4double y = boxYMin2 + (i / 48) * yStep2;
    G4double z = 0.5 * (FM1_5High + BM1_5High);
    G4ThreeVector posTSV_BBVIA30_FM1_5_BM1_5(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA30_FM1_5_BM1_5,
                      logicTSV_BBVIA30_FM1_5_BM1_5,
                      "TSV_BBVIA30_FM1_5_BM1_5",
                      logicEnv,
                      false,
                      i + 192,
                      true);
}
G4int numTSVs3 = 148;
G4double boxXMin3 = -5.285 * mm + TSV_BBVIA30_FM1_5_BM1_5R;
G4double boxXMax3 = -4.5 * mm - TSV_BBVIA30_FM1_5_BM1_5R;
G4double boxYMin3 = -4.065 * mm + 5*TSV_BBVIA30_FM1_5_BM1_5R;
G4double boxYMax3 = 4.065 * mm - 5*TSV_BBVIA30_FM1_5_BM1_5R;
G4double xStep3 = (boxXMax3 - boxXMin3) / 3;//4行
G4double yStep3 = (boxYMax3 - boxYMin3) / 36;//37列
for (int i = 0; i < numTSVs3; ++i) {
    G4double x = boxXMin3 + (i % 4) * xStep3;
    G4double y = boxYMin3 + (i / 4) * yStep3;
    G4double z = 0.5 * (FM1_5High + BM1_5High);
    G4ThreeVector posTSV_BBVIA30_FM1_5_BM1_5(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA30_FM1_5_BM1_5,
                      logicTSV_BBVIA30_FM1_5_BM1_5,
                      "TSV_BBVIA30_FM1_5_BM1_5",
                      logicEnv,
                      false,
                      i + 384,
                      true);
}
G4int numTSVs4 = 147;
G4double boxXMin4 = 4.5 * mm + TSV_BBVIA30_FM1_5_BM1_5R;
G4double boxXMax4 = 5.285 * mm - TSV_BBVIA30_FM1_5_BM1_5R;
G4double boxYMin4 = -4.065 * mm + 5*TSV_BBVIA30_FM1_5_BM1_5R;
G4double boxYMax4 = 4.065 * mm - 5*TSV_BBVIA30_FM1_5_BM1_5R;
G4double xStep4 = (boxXMax4 - boxXMin4) / 3;
G4double yStep4 = (boxYMax4 - boxYMin4) / 36;
for (int i = 0; i < numTSVs4; ++i) {
    G4double x = boxXMin4 + (i % 4) * xStep4;
    G4double y = boxYMin4 + (i / 4) * yStep4;
    G4double z = 0.5 * (FM1_5High + BM1_5High);
    G4ThreeVector posTSV_BBVIA30_FM1_5_BM1_5(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA30_FM1_5_BM1_5,
                      logicTSV_BBVIA30_FM1_5_BM1_5,
                      "TSV_BBVIA30_FM1_5_BM1_5",
                      logicEnv,
                      false,
                      i + 532,
                      true);
}

//放置679个TSV_BBVIA30_FM1_5_BM1_5,abs(x)>4.5*mm,abs(y)>4.065*mm

// for (int i = 0; i < 679; ++i) {
//         G4double x = (G4UniformRand() - 0.5) * (10.57*mm - 2 * TSV_BBVIA30_FM1_5_BM1_5R);
//         G4double y = (G4UniformRand() - 0.5) * (9.7*mm - 2 * TSV_BBVIA30_FM1_5_BM1_5R);
//     while(abs(x) < 4.5*mm +TSV_BBVIA30_FM1_5_BM1_5R && abs(y) < 4.065*mm +TSV_BBVIA30_FM1_5_BM1_5R){
//         G4double x = (G4UniformRand() - 0.5) * (10.57*mm - 2 * TSV_BBVIA30_FM1_5_BM1_5R);
//         G4double y = (G4UniformRand() - 0.5) * (9.7*mm - 2 * TSV_BBVIA30_FM1_5_BM1_5R);
//     }
//     G4double z = 0.5 * (BM1_5High + BM2_5High);
//     G4ThreeVector posTSV_BBVIA30_FM1_5_BM1_5(x, y, z);
//     while (OverlapWithOtherTSV_BBVIA30_FM1_5_BM1_5(posTSV_BBVIA30_FM1_5_BM1_5, TSV_BBVIA30_FM1_5_BM1_5R)) {
//         G4double x = (G4UniformRand() - 0.5) * (10.57*mm - 2 * TSV_BBVIA30_FM1_5_BM1_5R);
//         G4double y = (G4UniformRand() - 0.5) * (9.7*mm - 2 * TSV_BBVIA30_FM1_5_BM1_5R);
//         while(abs(x) < 4.5*mm +TSV_BBVIA30_FM1_5_BM1_5R && abs(y) < 4.065*mm +TSV_BBVIA30_FM1_5_BM1_5R){
//             G4double x = (G4UniformRand() - 0.5) * (10.57*mm - 2 * TSV_BBVIA30_FM1_5_BM1_5R);
//             G4double y = (G4UniformRand() - 0.5) * (9.7*mm - 2 * TSV_BBVIA30_FM1_5_BM1_5R);
//         }
//         posTSV_BBVIA30_FM1_5_BM1_5 = G4ThreeVector(x, y, z);
//     }
//     new G4PVPlacement(0,
//                       posTSV_BBVIA30_FM1_5_BM1_5,
//                       logicTSV_BBVIA30_FM1_5_BM1_5,
//                       "TSV_BBVIA30_FM1_5_BM1_5",
//                       logicEnv,
//                       false,
//                       i,
//                       true);
// }



//FP1_5,type:Dielectric,Material:FR_4
G4double FP1_5High = 0.42768*mm;
G4Box* solidFP1_5 =
    new G4Box("FP1_5",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FP1_5High-FM1_5High));
G4LogicalVolume* logicFP1_5 =
    new G4LogicalVolume(solidFP1_5,
                        EpoxyResin,
                        "FP1_5");
G4ThreeVector posFP1_5 = G4ThreeVector(0,0,0.5*(FP1_5High+FM1_5High));
    new G4PVPlacement(0,
                      posFP1_5,
                      logicFP1_5,
                      "FP1_5",
                      logicEnv,
                      false,
                      0,
                      true);                            

//UBM_5层,type:Metal,Material:Filled
G4double UBM_5High = 0.43168*mm;
G4Box* solidUBM_5 =
    new G4Box("UBM_5",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UBM_5High-FP1_5High));
G4LogicalVolume* logicUBM_5 =
    new G4LogicalVolume(solidUBM_5,
                        Filled,
                        "UBM_5");
G4ThreeVector posUBM_5 = G4ThreeVector(0,0,0.5*(UBM_5High+FP1_5High));
    new G4PVPlacement(0,
                      posUBM_5,
                      logicUBM_5,
                      "UBM_5",
                      logicEnv,
                      false,
                      0,
                      true);

//UNNAMED_054,type:Dielectric,Material:Air
G4double UNNAMED_054High = 0.51168*mm;
G4Box* solidUNNAMED_054 =
    new G4Box("UNNAMED_054",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UNNAMED_054High-UBM_5High));
G4LogicalVolume* logicUNNAMED_054 =
    new G4LogicalVolume(solidUNNAMED_054,
                        air,
                        "UNNAMED_054");
G4ThreeVector posUNNAMED_054 = G4ThreeVector(0,0,0.5*(UNNAMED_054High+UBM_5High));
    new G4PVPlacement(0,
                      posUNNAMED_054,
                      logicUNNAMED_054,
                      "UNNAMED_054",
                      logicEnv,
                      false,
                      0,
                      true);

//UBM_B_4,type:Metal,Material:Filled
G4double UBM_B_4High = 0.51568*mm;
G4Box* solidUBM_B_4 =
    new G4Box("UBM_B_4",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UBM_B_4High-UNNAMED_054High));
G4LogicalVolume* logicUBM_B_4 =
    new G4LogicalVolume(solidUBM_B_4,
                        Filled,
                        "UBM_B_4");
G4ThreeVector posUBM_B_4 = G4ThreeVector(0,0,0.5*(UBM_B_4High+UNNAMED_054High));
    new G4PVPlacement(0,
                      posUBM_B_4,
                      logicUBM_B_4,
                      "UBM_B_4",
                      logicEnv,
                      false,
                      0,
                      true);

//BP2_4,type:Dielectric,Material:FR_4
G4double BP2_4High = 0.71888*mm;
G4Box* solidBP2_4 =
    new G4Box("BP2_4",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BP2_4High-UBM_B_4High));
G4LogicalVolume* logicBP2_4 =
    new G4LogicalVolume(solidBP2_4,
                        EpoxyResin,
                        "BP2_4");
G4ThreeVector posBP2_4 = G4ThreeVector(0,0,0.5*(BP2_4High+UBM_B_4High));
    new G4PVPlacement(0,
                      posBP2_4,
                      logicBP2_4,
                      "BP2_4",
                      logicEnv,
                      false,
                      0,
                      true);  

//BM2_4层,type:Metal,Material:Filled
G4double BM2_4High = 0.74936*mm;
G4Box* solidBM2_4 =
    new G4Box("BM2_4",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BM2_4High-BP2_4High));
G4LogicalVolume* logicBM2_4 =
    new G4LogicalVolume(solidBM2_4,
                        Filled,
                        "BM2_4");
G4ThreeVector posBM2_4 = G4ThreeVector(0,0,0.5*(BM2_4High+BP2_4High));
    new G4PVPlacement(0,
                      posBM2_4,
                      logicBM2_4,
                      "BM2_4",
                      logicEnv,
                      false,
                      0,
                      true);

//焊盘：TSVINNERBGAPAD4,BM2_4-FM1_5,r=0.03mm
G4double TSVINNERBGAPAD4R = 0.03*mm;
G4Tubs* solidTSVINNERBGAPAD4 =
    new G4Tubs("TSVINNERBGAPAD4",
    0,
    TSVINNERBGAPAD4R,
    0.5*(BM2_4High-FM1_5High),
    0.*deg,
    360.*deg);
G4LogicalVolume* logicTSVINNERBGAPAD4 =
    new G4LogicalVolume(solidTSVINNERBGAPAD4,
                        copper,
                        "TSVINNERBGAPAD4");
G4int nTSVINNERBGAPAD4 = 0;                        
for(G4double x = -4*mm; x<=4*mm; x+=0.5*mm){
    for(G4double y = 0.625*mm; y<=3.625*mm; y+=0.5*mm){
        G4double z = 0.5 * (FM1_5High + BM2_4High);
        G4ThreeVector posTSVINNERBGAPAD4 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD4,
                      logicTSVINNERBGAPAD4,
                      "TSVINNERBGAPAD4",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD4,
                      true);
        nTSVINNERBGAPAD4++;
    }
}
for(G4double x = -4*mm; x<=4*mm; x+=0.5*mm){
    for(G4double y = -3.625*mm; y<=-0.125*mm; y+=0.5*mm){
        G4double z = 0.5 * (FM1_5High + BM2_4High);
        G4ThreeVector posTSVINNERBGAPAD4 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD4,
                      logicTSVINNERBGAPAD4,
                      "TSVINNERBGAPAD4",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD4,
                      true);
        nTSVINNERBGAPAD4++;
    }
}
for(G4double x = -5.125*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = 4.125*mm; y<=4.625*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_5High + BM2_4High);
        G4ThreeVector posTSVINNERBGAPAD4 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD4,
                      logicTSVINNERBGAPAD4,
                      "TSVINNERBGAPAD4",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD4,
                      true);
        nTSVINNERBGAPAD4++;
    }
}

for(G4double x = -5.125*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = -4.625*mm; y<=-4.125*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_5High + BM2_4High);
        G4ThreeVector posTSVINNERBGAPAD4 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD4,
                      logicTSVINNERBGAPAD4,
                      "TSVINNERBGAPAD4",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD4,
                      true);
        nTSVINNERBGAPAD4++;
    }
}

for(G4double x = -5.125*mm; x<=-4.625*mm; x+=0.25*mm){
    for(G4double y = -3.875*mm; y<=3.875*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_5High + BM2_4High);
        G4ThreeVector posTSVINNERBGAPAD4 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD4,
                      logicTSVINNERBGAPAD4,
                      "TSVINNERBGAPAD4",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD4,
                      true);
        nTSVINNERBGAPAD4++;
    }
}

for(G4double x = 4.625*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = -3.875*mm; y<=3.875*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_5High + BM2_4High);
        G4ThreeVector posTSVINNERBGAPAD4 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD4,
                      logicTSVINNERBGAPAD4,
                      "TSVINNERBGAPAD4",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD4,
                      true);
        nTSVINNERBGAPAD4++;
    }
}



//BP1_4层,type:Dielectric,Material:FR_4
G4double BP1_4High = 0.75736*mm;
G4Box* solidBP1_4 =
    new G4Box("BP1_4",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BP1_4High-BM2_4High));
G4LogicalVolume* logicBP1_4 =
    new G4LogicalVolume(solidBP1_4,
                        EpoxyResin,
                        "BP1_4");
G4ThreeVector posBP1_4 = G4ThreeVector(0,0,0.5*(BP1_4High+BM2_4High));
    new G4PVPlacement(0,
                      posBP1_4,
                      logicBP1_4,
                      "BP1_4",
                      logicEnv,
                      false,
                      0,
                      true);

//BM1_4层,type:Metal,Material:Filled
G4double BM1_4High = 0.76136*mm;
G4Box* solidBM1_4 =
    new G4Box("BM1_4",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BM1_4High-BP1_4High));
G4LogicalVolume* logicBM1_4 =
    new G4LogicalVolume(solidBM1_4,
                        Filled,
                        "BM1_4");
G4ThreeVector posBM1_4 = G4ThreeVector(0,0,0.5*(BM1_4High+BP1_4High));
    new G4PVPlacement(0,
                      posBM1_4,
                      logicBM1_4,
                      "BM1_4",
                      logicEnv,
                      false,
                      0,
                      true);


//焊盘：TSV-BBVIA30-BM1_4_BM2_4,copper,897，0.015mm
G4double TSV_BBVIA30_BM1_4_BM2_4R = 0.015*mm;
G4Tubs* solidTSV_BBVIA30_BM1_4_BM2_4 =
    new G4Tubs("TSV_BBVIA30_BM1_4_BM2_4",
    0,
    TSV_BBVIA30_BM1_4_BM2_4R,
    0.5*(BM1_4High - BM2_4High),
    0.*deg,
    360.*deg);
G4LogicalVolume* logicTSV_BBVIA30_BM1_4_BM2_4 =
    new G4LogicalVolume(solidTSV_BBVIA30_BM1_4_BM2_4,
                        copper,
                        "TSV_BBVIA30_BM1_4_BM2_4");
G4double TSV_BBVIA30_BM1_4_BM2_4Xmin1 = 0.2*mm + TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Xmax1 = 4.3*mm - TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Ymin1 = 0.4*mm + TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Ymax1 = 3.9*mm - TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Xstep1 = (TSV_BBVIA30_BM1_4_BM2_4Xmax1 - TSV_BBVIA30_BM1_4_BM2_4Xmin1)/15;
G4double TSV_BBVIA30_BM1_4_BM2_4Ystep1 = (TSV_BBVIA30_BM1_4_BM2_4Ymax1 - TSV_BBVIA30_BM1_4_BM2_4Ymin1)/11;
for(int nTSV_BBVIA30_BM1_4_BM2_4=0; nTSV_BBVIA30_BM1_4_BM2_4<132; nTSV_BBVIA30_BM1_4_BM2_4++){
    G4double z = 0.5 * (BM1_4High + BM2_4High);
    G4double x = TSV_BBVIA30_BM1_4_BM2_4Xmin1 + TSV_BBVIA30_BM1_4_BM2_4Xstep1*(nTSV_BBVIA30_BM1_4_BM2_4%16);
    G4double y = TSV_BBVIA30_BM1_4_BM2_4Ymin1 + TSV_BBVIA30_BM1_4_BM2_4Ystep1*(nTSV_BBVIA30_BM1_4_BM2_4/16);
    G4ThreeVector posTSV_BBVIA30_BM1_4_BM2_4 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_4_BM2_4,
                      logicTSV_BBVIA30_BM1_4_BM2_4,
                      "TSV_BBVIA30_BM1_4_BM2_4",
                      logicEnv,
                      false,
                      nTSV_BBVIA30_BM1_4_BM2_4,
                      true);
}
G4double TSV_BBVIA30_BM1_4_BM2_4Xmin2 = -4.2*mm + TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Xmax2 = 4.2*mm - TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Ymin2 = -3.8*mm + TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Ymax2 = - TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Xstep2 = (TSV_BBVIA30_BM1_4_BM2_4Xmax2 - TSV_BBVIA30_BM1_4_BM2_4Xmin2)/36;
G4double TSV_BBVIA30_BM1_4_BM2_4Ystep2 = (TSV_BBVIA30_BM1_4_BM2_4Ymax2 - TSV_BBVIA30_BM1_4_BM2_4Ymin2)/3;
for(int nTSV_BBVIA30_BM1_4_BM2_4=0; nTSV_BBVIA30_BM1_4_BM2_4<148; nTSV_BBVIA30_BM1_4_BM2_4++){
    G4double z = 0.5 * (BM1_4High + BM2_4High);
    G4double x = TSV_BBVIA30_BM1_4_BM2_4Xmin2 + TSV_BBVIA30_BM1_4_BM2_4Xstep2*(nTSV_BBVIA30_BM1_4_BM2_4%37);
    G4double y = TSV_BBVIA30_BM1_4_BM2_4Ymin2 + TSV_BBVIA30_BM1_4_BM2_4Ystep2*(nTSV_BBVIA30_BM1_4_BM2_4/37);
    G4ThreeVector posTSV_BBVIA30_BM1_4_BM2_4 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_4_BM2_4,
                      logicTSV_BBVIA30_BM1_4_BM2_4,
                      "TSV_BBVIA30_BM1_4_BM2_4",
                      logicEnv,
                      false,
                      nTSV_BBVIA30_BM1_4_BM2_4+132,
                      true);
}
G4double TSV_BBVIA30_BM1_4_BM2_4Xmin3 = -5.3*mm + TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Xmax3 = 5.3*mm - TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Ymin3 = -4.8*mm + TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Ymax3 = -4.0*mm - 5*TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Xstep3 = (TSV_BBVIA30_BM1_4_BM2_4Xmax3 - TSV_BBVIA30_BM1_4_BM2_4Xmin3)/43;
G4double TSV_BBVIA30_BM1_4_BM2_4Ystep3 = (TSV_BBVIA30_BM1_4_BM2_4Ymax3 - TSV_BBVIA30_BM1_4_BM2_4Ymin3)/3;
for (int nTSV_BBVIA30_BM1_4_BM2_4=0; nTSV_BBVIA30_BM1_4_BM2_4<176; nTSV_BBVIA30_BM1_4_BM2_4++){
    G4double z = 0.5 * (BM1_4High + BM2_4High);
    G4double x = TSV_BBVIA30_BM1_4_BM2_4Xmin3 + TSV_BBVIA30_BM1_4_BM2_4Xstep3*(nTSV_BBVIA30_BM1_4_BM2_4%44);
    G4double y = TSV_BBVIA30_BM1_4_BM2_4Ymin3 + TSV_BBVIA30_BM1_4_BM2_4Ystep3*(nTSV_BBVIA30_BM1_4_BM2_4/44);
    G4ThreeVector posTSV_BBVIA30_BM1_4_BM2_4 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_4_BM2_4,
                      logicTSV_BBVIA30_BM1_4_BM2_4,
                      "TSV_BBVIA30_BM1_4_BM2_4",
                      logicEnv,
                      false,
                      nTSV_BBVIA30_BM1_4_BM2_4+280,
                      true);
}
G4double TSV_BBVIA30_BM1_4_BM2_4Xmin4 = -5.3*mm + TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Xmax4 = 5.3*mm - TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Ymin4 = 4.0*mm + 5*TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Ymax4 = 4.8*mm - TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Xstep4 = (TSV_BBVIA30_BM1_4_BM2_4Xmax4 - TSV_BBVIA30_BM1_4_BM2_4Xmin4)/43;
G4double TSV_BBVIA30_BM1_4_BM2_4Ystep4 = (TSV_BBVIA30_BM1_4_BM2_4Ymax4 - TSV_BBVIA30_BM1_4_BM2_4Ymin4)/3;
for (int nTSV_BBVIA30_BM1_4_BM2_4=0; nTSV_BBVIA30_BM1_4_BM2_4<176; nTSV_BBVIA30_BM1_4_BM2_4++){
    G4double z = 0.5 * (BM1_4High + BM2_4High);
    G4double x = TSV_BBVIA30_BM1_4_BM2_4Xmin4 + TSV_BBVIA30_BM1_4_BM2_4Xstep4*(nTSV_BBVIA30_BM1_4_BM2_4%44);
    G4double y = TSV_BBVIA30_BM1_4_BM2_4Ymin4 + TSV_BBVIA30_BM1_4_BM2_4Ystep4*(nTSV_BBVIA30_BM1_4_BM2_4/44);
    G4ThreeVector posTSV_BBVIA30_BM1_4_BM2_4 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_4_BM2_4,
                      logicTSV_BBVIA30_BM1_4_BM2_4,
                      "TSV_BBVIA30_BM1_4_BM2_4",
                      logicEnv,
                      false,
                      nTSV_BBVIA30_BM1_4_BM2_4+456,
                      true);
}
G4double TSV_BBVIA30_BM1_4_BM2_4Xmin5 = -5.3*mm + TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Xmax5 = -4.5*mm - TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Ymin5 = -4.0*mm + 5*TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Ymax5 = 4.0*mm - 5*TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Xstep5 = (TSV_BBVIA30_BM1_4_BM2_4Xmax5 - TSV_BBVIA30_BM1_4_BM2_4Xmin5)/4;
G4double TSV_BBVIA30_BM1_4_BM2_4Ystep5 = (TSV_BBVIA30_BM1_4_BM2_4Ymax5 - TSV_BBVIA30_BM1_4_BM2_4Ymin5)/26;
for (int nTSV_BBVIA30_BM1_4_BM2_4=0; nTSV_BBVIA30_BM1_4_BM2_4<133; nTSV_BBVIA30_BM1_4_BM2_4++){
    G4double z = 0.5 * (BM1_4High + BM2_4High);
    G4double x = TSV_BBVIA30_BM1_4_BM2_4Xmin5 + TSV_BBVIA30_BM1_4_BM2_4Xstep5*(nTSV_BBVIA30_BM1_4_BM2_4%5);
    G4double y = TSV_BBVIA30_BM1_4_BM2_4Ymin5 + TSV_BBVIA30_BM1_4_BM2_4Ystep5*(nTSV_BBVIA30_BM1_4_BM2_4/5);
    G4ThreeVector posTSV_BBVIA30_BM1_4_BM2_4 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_4_BM2_4,
                      logicTSV_BBVIA30_BM1_4_BM2_4,
                      "TSV_BBVIA30_BM1_4_BM2_4",
                      logicEnv,
                      false,
                      nTSV_BBVIA30_BM1_4_BM2_4+632,
                      true);
}
G4double TSV_BBVIA30_BM1_4_BM2_4Xmin6 = 4.5*mm + TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Xmax6 = 5.3*mm - TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Ymin6 = -4.0*mm + 5*TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Ymax6 = 4.0*mm - 5*TSV_BBVIA30_BM1_4_BM2_4R;
G4double TSV_BBVIA30_BM1_4_BM2_4Xstep6 = (TSV_BBVIA30_BM1_4_BM2_4Xmax6 - TSV_BBVIA30_BM1_4_BM2_4Xmin6)/4;
G4double TSV_BBVIA30_BM1_4_BM2_4Ystep6 = (TSV_BBVIA30_BM1_4_BM2_4Ymax6 - TSV_BBVIA30_BM1_4_BM2_4Ymin6)/26;
for (int nTSV_BBVIA30_BM1_4_BM2_4=0; nTSV_BBVIA30_BM1_4_BM2_4<132; nTSV_BBVIA30_BM1_4_BM2_4++){
    G4double z = 0.5 * (BM1_4High + BM2_4High);
    G4double x = TSV_BBVIA30_BM1_4_BM2_4Xmin6 + TSV_BBVIA30_BM1_4_BM2_4Xstep6*(nTSV_BBVIA30_BM1_4_BM2_4%5);
    G4double y = TSV_BBVIA30_BM1_4_BM2_4Ymin6 + TSV_BBVIA30_BM1_4_BM2_4Ystep6*(nTSV_BBVIA30_BM1_4_BM2_4/5);
    G4ThreeVector posTSV_BBVIA30_BM1_4_BM2_4 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_4_BM2_4,
                      logicTSV_BBVIA30_BM1_4_BM2_4,
                      "TSV_BBVIA30_BM1_4_BM2_4",
                      logicEnv,
                      false,
                      nTSV_BBVIA30_BM1_4_BM2_4+765,
                      true);
}





//BPI_CORE_4,type:Dielectric,Material:FR_4
G4double BPI_CORE_4High = 0.76936*mm;
G4Box* solidBPI_CORE_4 =
    new G4Box("BPI_CORE_4",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BPI_CORE_4High-BM1_4High));
G4LogicalVolume* logicBPI_CORE_4 =
    new G4LogicalVolume(solidBPI_CORE_4,
                        EpoxyResin,
                        "BPI_CORE_4");
G4ThreeVector posBPI_CORE_4 = G4ThreeVector(0,0,0.5*(BPI_CORE_4High+BM1_4High));
    new G4PVPlacement(0,
                      posBPI_CORE_4,
                      logicBPI_CORE_4,
                      "BPI_CORE_4",
                      logicEnv,
                      false,
                      0,
                      true);

//CORE2_4,type:Dielectric,Material:FR_4
G4double CORE2_4High = 0.79936*mm;
G4Box* solidCORE2_4 =
    new G4Box("CORE2_4",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(CORE2_4High-BPI_CORE_4High));
G4LogicalVolume* logicCORE2_4 =
    new G4LogicalVolume(solidCORE2_4,
                        EpoxyResin,
                        "CORE2_4");
G4ThreeVector posCORE2_4 = G4ThreeVector(0,0,0.5*(CORE2_4High+BPI_CORE_4High));
    new G4PVPlacement(0,
                      posCORE2_4,
                      logicCORE2_4,
                      "CORE2_4",
                      logicEnv,
                      false,
                      0,
                      true);

//CORE1_4,type:Dielectric,Material:FR_4
G4double CORE1_4High = 0.91936*mm;      
G4Box* solidCORE1_4 =
    new G4Box("CORE1_4",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(CORE1_4High-CORE2_4High));
G4LogicalVolume* logicCORE1_4 =
    new G4LogicalVolume(solidCORE1_4,
                        EpoxyResin,
                        "CORE1_4");
G4ThreeVector posCORE1_4 = G4ThreeVector(0,0,0.5*(CORE1_4High+CORE2_4High));
    new G4PVPlacement(0,
                      posCORE1_4,
                      logicCORE1_4,
                      "CORE1_4",
                      logicEnv,
                      false,
                      0,
                      true);

//FPI_CORE_4,type:Dielectric,Material:FR_4
G4double FPI_CORE_4High = 0.92736*mm;
G4Box* solidFPI_CORE_4 =
    new G4Box("FPI_CORE_4",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FPI_CORE_4High-CORE1_4High));
G4LogicalVolume* logicFPI_CORE_4 =
    new G4LogicalVolume(solidFPI_CORE_4,
                        EpoxyResin,
                        "FPI_CORE_4");
G4ThreeVector posFPI_CORE_4 = G4ThreeVector(0,0,0.5*(FPI_CORE_4High+CORE1_4High));
    new G4PVPlacement(0,
                      posFPI_CORE_4,
                      logicFPI_CORE_4,
                      "FPI_CORE_4",
                      logicEnv,
                      false,
                      0,
                      true);

//FM1_4,type:Metal,Material:Filled
G4double FM1_4High = 0.93136*mm;    
G4Box* solidFM1_4 =
    new G4Box("FM1_4",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FM1_4High-FPI_CORE_4High));
G4LogicalVolume* logicFM1_4 =
    new G4LogicalVolume(solidFM1_4,
                        Filled,
                        "FM1_4");
G4ThreeVector posFM1_4 = G4ThreeVector(0,0,0.5*(FM1_4High+FPI_CORE_4High));
    new G4PVPlacement(0,
                      posFM1_4,
                      logicFM1_4,
                      "FM1_4",
                      logicEnv,
                      false,
                      0,
                      true);

//焊盘：TSV-BBVIA-30-FM1_4_BM1_4,0.01mm
G4double TSV_BBVIA_30_FM1_4_BM1_4R = 0.01*mm;
G4Tubs* solidTSV_BBVIA_30_FM1_4_BM1_4 =
    new G4Tubs("TSV_BBVIA_30_FM1_4_BM1_4",
               0,
               TSV_BBVIA_30_FM1_4_BM1_4R,
               0.5*(FM1_4High-BM1_4High),
               0.*deg,
               360.*deg);
G4LogicalVolume* logicTSV_BBVIA_30_FM1_4_BM1_4 =
    new G4LogicalVolume(solidTSV_BBVIA_30_FM1_4_BM1_4,
                        copper,
                        "TSV_BBVIA_30_FM1_4_BM1_4");
G4double TSV_BBVIA_30_FM1_4_BM1_4Xmin1 = -5.25*mm + TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Xmax1 = 5.25*mm - TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Ymin1 = -4.8*mm + TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Ymax1 = -4.05*mm - 5*TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Xstep1 = (TSV_BBVIA_30_FM1_4_BM1_4Xmax1 - TSV_BBVIA_30_FM1_4_BM1_4Xmin1)/47;
G4double TSV_BBVIA_30_FM1_4_BM1_4Ystep1 = (TSV_BBVIA_30_FM1_4_BM1_4Ymax1 - TSV_BBVIA_30_FM1_4_BM1_4Ymin1)/3;
for (int nTSV_BBVIA_30_FM1_4_BM1_4=0; nTSV_BBVIA_30_FM1_4_BM1_4<192; nTSV_BBVIA_30_FM1_4_BM1_4++){
    G4double z = 0.5 * (FM1_4High + BM1_4High);
    G4double x = TSV_BBVIA_30_FM1_4_BM1_4Xmin1 + TSV_BBVIA_30_FM1_4_BM1_4Xstep1*(nTSV_BBVIA_30_FM1_4_BM1_4%48);
    G4double y = TSV_BBVIA_30_FM1_4_BM1_4Ymin1 + TSV_BBVIA_30_FM1_4_BM1_4Ystep1*(nTSV_BBVIA_30_FM1_4_BM1_4/48);
    G4ThreeVector posTSV_BBVIA_30_FM1_4_BM1_4 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA_30_FM1_4_BM1_4,
                      logicTSV_BBVIA_30_FM1_4_BM1_4,
                      "TSV_BBVIA_30_FM1_4_BM1_4",
                      logicEnv,
                      false,
                      nTSV_BBVIA_30_FM1_4_BM1_4,
                      true);
}
G4double TSV_BBVIA_30_FM1_4_BM1_4Xmin2 = -5.25*mm + TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Xmax2 = 5.25*mm - TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Ymin2 = 4.05*mm + 5*TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Ymax2 = 4.8*mm - TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Xstep2 = (TSV_BBVIA_30_FM1_4_BM1_4Xmax2 - TSV_BBVIA_30_FM1_4_BM1_4Xmin2)/47;
G4double TSV_BBVIA_30_FM1_4_BM1_4Ystep2 = (TSV_BBVIA_30_FM1_4_BM1_4Ymax2 - TSV_BBVIA_30_FM1_4_BM1_4Ymin2)/3;
for (int nTSV_BBVIA_30_FM1_4_BM1_4=0; nTSV_BBVIA_30_FM1_4_BM1_4<191; nTSV_BBVIA_30_FM1_4_BM1_4++){
    G4double z = 0.5 * (FM1_4High + BM1_4High);
    G4double x = TSV_BBVIA_30_FM1_4_BM1_4Xmin2 + TSV_BBVIA_30_FM1_4_BM1_4Xstep2*(nTSV_BBVIA_30_FM1_4_BM1_4%48);
    G4double y = TSV_BBVIA_30_FM1_4_BM1_4Ymin2 + TSV_BBVIA_30_FM1_4_BM1_4Ystep2*(nTSV_BBVIA_30_FM1_4_BM1_4/48);
    G4ThreeVector posTSV_BBVIA_30_FM1_4_BM1_4 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA_30_FM1_4_BM1_4,
                      logicTSV_BBVIA_30_FM1_4_BM1_4,
                      "TSV_BBVIA_30_FM1_4_BM1_4",
                      logicEnv,
                      false,
                      nTSV_BBVIA_30_FM1_4_BM1_4+192,
                      true);
}
G4double TSV_BBVIA_30_FM1_4_BM1_4Xmin3 = -5.25*mm + TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Xmax3 = -4.5*mm - TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Ymin3 = -4.05*mm + 5*TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Ymax3 = 4.05*mm - 5*TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Xstep3 = (TSV_BBVIA_30_FM1_4_BM1_4Xmax3 - TSV_BBVIA_30_FM1_4_BM1_4Xmin3)/3;
G4double TSV_BBVIA_30_FM1_4_BM1_4Ystep3 = (TSV_BBVIA_30_FM1_4_BM1_4Ymax3 - TSV_BBVIA_30_FM1_4_BM1_4Ymin3)/37;
for (int nTSV_BBVIA_30_FM1_4_BM1_4=0; nTSV_BBVIA_30_FM1_4_BM1_4<152; nTSV_BBVIA_30_FM1_4_BM1_4++){
    G4double z = 0.5 * (FM1_4High + BM1_4High);
    G4double x = TSV_BBVIA_30_FM1_4_BM1_4Xmin3 + TSV_BBVIA_30_FM1_4_BM1_4Xstep3*(nTSV_BBVIA_30_FM1_4_BM1_4%4);
    G4double y = TSV_BBVIA_30_FM1_4_BM1_4Ymin3 + TSV_BBVIA_30_FM1_4_BM1_4Ystep3*(nTSV_BBVIA_30_FM1_4_BM1_4/4);
    G4ThreeVector posTSV_BBVIA_30_FM1_4_BM1_4 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA_30_FM1_4_BM1_4,
                      logicTSV_BBVIA_30_FM1_4_BM1_4,
                      "TSV_BBVIA_30_FM1_4_BM1_4",
                      logicEnv,
                      false,
                      nTSV_BBVIA_30_FM1_4_BM1_4+383,
                      true);
}
G4double TSV_BBVIA_30_FM1_4_BM1_4Xmin4 = 4.5*mm + TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Xmax4 = 5.25*mm - TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Ymin4 = -4.05*mm + 5*TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Ymax4 = 4.05*mm - 5*TSV_BBVIA_30_FM1_4_BM1_4R;
G4double TSV_BBVIA_30_FM1_4_BM1_4Xstep4 = (TSV_BBVIA_30_FM1_4_BM1_4Xmax4 - TSV_BBVIA_30_FM1_4_BM1_4Xmin4)/3;
G4double TSV_BBVIA_30_FM1_4_BM1_4Ystep4 = (TSV_BBVIA_30_FM1_4_BM1_4Ymax4 - TSV_BBVIA_30_FM1_4_BM1_4Ymin4)/37;
for (int nTSV_BBVIA_30_FM1_4_BM1_4=0; nTSV_BBVIA_30_FM1_4_BM1_4<152; nTSV_BBVIA_30_FM1_4_BM1_4++){
    G4double z = 0.5 * (FM1_4High + BM1_4High);
    G4double x = TSV_BBVIA_30_FM1_4_BM1_4Xmin4 + TSV_BBVIA_30_FM1_4_BM1_4Xstep4*(nTSV_BBVIA_30_FM1_4_BM1_4%4);
    G4double y = TSV_BBVIA_30_FM1_4_BM1_4Ymin4 + TSV_BBVIA_30_FM1_4_BM1_4Ystep4*(nTSV_BBVIA_30_FM1_4_BM1_4/4);
    G4ThreeVector posTSV_BBVIA_30_FM1_4_BM1_4 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA_30_FM1_4_BM1_4,
                      logicTSV_BBVIA_30_FM1_4_BM1_4,
                      "TSV_BBVIA_30_FM1_4_BM1_4",
                      logicEnv,
                      false,
                      nTSV_BBVIA_30_FM1_4_BM1_4+535,
                      true);
}



//FP1_4,type:Dielectric,Material:FR_4
G4double FP1_4High = 0.93936*mm;
G4Box* solidFP1_4 =
    new G4Box("FP1_4",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FP1_4High-FM1_4High));
G4LogicalVolume* logicFP1_4 =
    new G4LogicalVolume(solidFP1_4,
                        EpoxyResin,
                        "FP1_4");
G4ThreeVector posFP1_4 = G4ThreeVector(0,0,0.5*(FP1_4High+FM1_4High));
    new G4PVPlacement(0,
                      posFP1_4,
                      logicFP1_4,
                      "FP1_4",
                      logicEnv,
                      false,
                      0,
                      true);

//UBM_4,type:Metal,Material:Filled
G4double UBM_4High = 0.94336*mm;
G4Box* solidUBM_4 =
    new G4Box("UBM_4",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UBM_4High-FP1_4High));
G4LogicalVolume* logicUBM_4 =
    new G4LogicalVolume(solidUBM_4,
                        Filled,
                        "UBM_4");
G4ThreeVector posUBM_4 = G4ThreeVector(0,0,0.5*(UBM_4High+FP1_4High));
    new G4PVPlacement(0,
                      posUBM_4,
                      logicUBM_4,
                      "UBM_4",
                      logicEnv,
                      false,
                      0,
                      true);

//UNNAMED_041,type:Dielectric,Material:Air
G4double UNNAMED_041High = 0.99336*mm;
G4Box* solidUNNAMED_041 =
    new G4Box("UNNAMED_041",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UNNAMED_041High-UBM_4High));
G4LogicalVolume* logicUNNAMED_041 =
    new G4LogicalVolume(solidUNNAMED_041,
                        air,
                        "UNNAMED_041");
G4ThreeVector posUNNAMED_041 = G4ThreeVector(0,0,0.5*(UNNAMED_041High+UBM_4High));
    new G4PVPlacement(0,
                      posUNNAMED_041,
                      logicUNNAMED_041,
                      "UNNAMED_041",
                      logicEnv,
                      false,
                      0,
                      true);

//UBM_B_3,type:Metal,Material:Filled
G4double UBM_B_3High = 0.99736*mm;
G4Box* solidUBM_B_3 =
    new G4Box("UBM_B_3",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UBM_B_3High-UNNAMED_041High));
G4LogicalVolume* logicUBM_B_3 =
    new G4LogicalVolume(solidUBM_B_3,
                        Filled,
                        "UBM_B_3");
G4ThreeVector posUBM_B_3 = G4ThreeVector(0,0,0.5*(UBM_B_3High+UNNAMED_041High));
    new G4PVPlacement(0,
                      posUBM_B_3,
                      logicUBM_B_3,
                      "UBM_B_3",
                      logicEnv,
                      false,
                      0,
                      true);

//BP2_3层,type:Dielectric,Material:FR_4
G4double BP2_3High = 1.20056*mm;
G4Box* solidBP2_3 =
    new G4Box("BP2_3",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BP2_3High-UBM_B_3High));
G4LogicalVolume* logicBP2_3 =
    new G4LogicalVolume(solidBP2_3,
                        EpoxyResin,
                        "BP2_3");
G4ThreeVector posBP2_3 = G4ThreeVector(0,0,0.5*(BP2_3High+UBM_B_3High));
    new G4PVPlacement(0,
                      posBP2_3,
                      logicBP2_3,
                      "BP2_3",
                      logicEnv,
                      false,
                      0,
                      true);

//BM2_3,type:Metal,Material:Filled
G4double BM2_3High = 1.23104*mm;
G4Box* solidBM2_3 =
    new G4Box("BM2_3",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BM2_3High-BP2_3High));
G4LogicalVolume* logicBM2_3 =
    new G4LogicalVolume(solidBM2_3,
                        Filled,
                        "BM2_3");
G4ThreeVector posBM2_3 = G4ThreeVector(0,0,0.5*(BM2_3High+BP2_3High));
    new G4PVPlacement(0,
                      posBM2_3,
                      logicBM2_3,
                      "BM2_3",
                      logicEnv,
                      false,
                      0,
                      true);

//焊盘：TSVINNERBGAPAD3,0.03*mm,copper
G4double TSVINNERBGAPAD3R = 0.03*mm;
G4Tubs* solidTSVINNERBGAPAD3 =
    new G4Tubs("TSVINNERBGAPAD3",
    0,
    TSVINNERBGAPAD3R,
    0.5*(BM2_3High-FM1_4High),
    0.*deg,
    360.*deg);
G4LogicalVolume* logicTSVINNERBGAPAD3 =
    new G4LogicalVolume(solidTSVINNERBGAPAD3,
                        copper,
                        "TSVINNERBGAPAD3");
G4int nTSVINNERBGAPAD3 = 0;                        
for(G4double x = -4*mm; x<=4*mm; x+=0.5*mm){
    for(G4double y = 0.625*mm; y<=3.625*mm; y+=0.5*mm){
        G4double z = 0.5 * (FM1_4High + BM2_3High);
        G4ThreeVector posTSVINNERBGAPAD3 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD3,
                      logicTSVINNERBGAPAD3,
                      "TSVINNERBGAPAD3",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD3,
                      true);
        nTSVINNERBGAPAD3++;
    }
}
for(G4double x = -4*mm; x<=4*mm; x+=0.5*mm){
    for(G4double y = -3.625*mm; y<=-0.125*mm; y+=0.5*mm){
        G4double z = 0.5 * (FM1_4High + BM2_3High);
        G4ThreeVector posTSVINNERBGAPAD3 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD3,
                      logicTSVINNERBGAPAD3,
                      "TSVINNERBGAPAD3",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD3,
                      true);
        nTSVINNERBGAPAD3++;
    }
}
for(G4double x = -5.125*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = 4.125*mm; y<=4.625*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_4High + BM2_3High);
        G4ThreeVector posTSVINNERBGAPAD3 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD3,
                      logicTSVINNERBGAPAD3,
                      "TSVINNERBGAPAD3",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD3,
                      true);
        nTSVINNERBGAPAD3++;
    }
}

for(G4double x = -5.125*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = -4.625*mm; y<=-4.125*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_4High + BM2_3High);
        G4ThreeVector posTSVINNERBGAPAD3 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD3,
                      logicTSVINNERBGAPAD3,
                      "TSVINNERBGAPAD3",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD3,
                      true);
        nTSVINNERBGAPAD3++;
    }
}

for(G4double x = -5.125*mm; x<=-4.625*mm; x+=0.25*mm){
    for(G4double y = -3.875*mm; y<=3.875*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_4High + BM2_3High);
        G4ThreeVector posTSVINNERBGAPAD3 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD3,
                      logicTSVINNERBGAPAD3,
                      "TSVINNERBGAPAD3",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD3,
                      true);
        nTSVINNERBGAPAD3++;
    }
}

for(G4double x = 4.625*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = -3.875*mm; y<=3.875*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_4High + BM2_3High);
        G4ThreeVector posTSVINNERBGAPAD3 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD3,
                      logicTSVINNERBGAPAD3,
                      "TSVINNERBGAPAD3",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD3,
                      true);
        nTSVINNERBGAPAD3++;
    }
}

//BP1_3,type:Dielectric,Material:FR_4
G4double BP1_3High = 1.23904*mm;
G4Box* solidBP1_3 =
    new G4Box("BP1_3",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BP1_3High-BM2_3High));
G4LogicalVolume* logicBP1_3 =
    new G4LogicalVolume(solidBP1_3,
                        EpoxyResin,
                        "BP1_3");
G4ThreeVector posBP1_3 = G4ThreeVector(0,0,0.5*(BP1_3High+BM2_3High));
    new G4PVPlacement(0,
                      posBP1_3,
                      logicBP1_3,
                      "BP1_3",
                      logicEnv,
                      false,
                      0,
                      true);
                      
//BM1_3,type:Metal,Material:Filled
G4double BM1_3High = 1.24304*mm;  
G4Box* solidBM1_3 =
    new G4Box("BM1_3",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BM1_3High-BP1_3High));
G4LogicalVolume* logicBM1_3 =
    new G4LogicalVolume(solidBM1_3,
                        Filled,
                        "BM1_3");
G4ThreeVector posBM1_3 = G4ThreeVector(0,0,0.5*(BM1_3High+BP1_3High));
    new G4PVPlacement(0,
                      posBM1_3,
                      logicBM1_3,
                      "BM1_3",
                      logicEnv,
                      false,
                      0,
                      true);
//焊盘，TSV-BBVIA30-BM1_3_BM2_3,0.0015mm,copper
G4double TSV_BBVIA30_BM1_3_BM2_3R = 0.015*mm;
G4Tubs* solidTSV_BBVIA30_BM1_3_BM2_3 =
    new G4Tubs("TSV_BBVIA30_BM1_3_BM2_3",
    0,
    TSV_BBVIA30_BM1_3_BM2_3R,
    0.5*(BM1_3High-BM2_3High),
    0.*deg,
    360.*deg);
G4LogicalVolume* logicTSV_BBVIA30_BM1_3_BM2_3 =
    new G4LogicalVolume(solidTSV_BBVIA30_BM1_3_BM2_3,
                        copper,
                        "TSV_BBVIA30_BM1_3_BM2_3");
G4double TSV_BBVIA30_BM1_3_BM2_3Xmin1 = 0.4*mm + TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Xmax1 = 3.6*mm - TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Ymin1 = 0.5*mm + TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Ymax1 = 3.75*mm - TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Xstep1 = (TSV_BBVIA30_BM1_3_BM2_3Xmax1 - TSV_BBVIA30_BM1_3_BM2_3Xmin1)/9;
G4double TSV_BBVIA30_BM1_3_BM2_3Ystep1 = (TSV_BBVIA30_BM1_3_BM2_3Ymax1 - TSV_BBVIA30_BM1_3_BM2_3Ymin1)/9;
for(G4int nTSV_BBVIA30_BM1_3_BM2_3X = 0; nTSV_BBVIA30_BM1_3_BM2_3X < 100; nTSV_BBVIA30_BM1_3_BM2_3X++){
    G4double x = TSV_BBVIA30_BM1_3_BM2_3Xmin1 + (nTSV_BBVIA30_BM1_3_BM2_3X%10)*TSV_BBVIA30_BM1_3_BM2_3Xstep1;
    G4double y = TSV_BBVIA30_BM1_3_BM2_3Ymin1 + (nTSV_BBVIA30_BM1_3_BM2_3X/10)*TSV_BBVIA30_BM1_3_BM2_3Ystep1;
    G4ThreeVector posTSV_BBVIA30_BM1_3_BM2_3 = G4ThreeVector(x, y, 0.5*(BM1_3High+BM2_3High));
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_3_BM2_3,
                      logicTSV_BBVIA30_BM1_3_BM2_3,
                      "TSV_BBVIA30_BM1_3_BM2_3",
                      logicEnv,
                      false,
                      nTSV_BBVIA30_BM1_3_BM2_3X,
                      true);                      
}
G4double TSV_BBVIA30_BM1_3_BM2_3step1 = 0.06126*mm;
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm + TSV_BBVIA30_BM1_3_BM2_3step1, 3.125*mm + TSV_BBVIA30_BM1_3_BM2_3step1, 0.5*(BM1_3High+BM2_3High)),
                logicTSV_BBVIA30_BM1_3_BM2_3,
                "TSV_BBVIA30_BM1_3_BM2_3",
                logicEnv,
                false,
                100,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm + TSV_BBVIA30_BM1_3_BM2_3step1, 3.125*mm - TSV_BBVIA30_BM1_3_BM2_3step1, 0.5*(BM1_3High+BM2_3High)),
                logicTSV_BBVIA30_BM1_3_BM2_3,
                "TSV_BBVIA30_BM1_3_BM2_3",
                logicEnv,
                false,
                101,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm - TSV_BBVIA30_BM1_3_BM2_3step1, 3.125*mm + TSV_BBVIA30_BM1_3_BM2_3step1, 0.5*(BM1_3High+BM2_3High)),
                logicTSV_BBVIA30_BM1_3_BM2_3,
                "TSV_BBVIA30_BM1_3_BM2_3",
                logicEnv,
                false,
                102,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm - TSV_BBVIA30_BM1_3_BM2_3step1, 3.125*mm - TSV_BBVIA30_BM1_3_BM2_3step1, 0.5*(BM1_3High+BM2_3High)),
                logicTSV_BBVIA30_BM1_3_BM2_3,
                "TSV_BBVIA30_BM1_3_BM2_3",
                logicEnv,
                false,
                103,
                true);    
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm + TSV_BBVIA30_BM1_3_BM2_3step1, 2.125*mm + TSV_BBVIA30_BM1_3_BM2_3step1, 0.5*(BM1_3High+BM2_3High)),
                logicTSV_BBVIA30_BM1_3_BM2_3,
                "TSV_BBVIA30_BM1_3_BM2_3",
                logicEnv,
                false,
                104,
                true); 
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm + TSV_BBVIA30_BM1_3_BM2_3step1, 2.125*mm - TSV_BBVIA30_BM1_3_BM2_3step1, 0.5*(BM1_3High+BM2_3High)),
                logicTSV_BBVIA30_BM1_3_BM2_3,
                "TSV_BBVIA30_BM1_3_BM2_3",
                logicEnv,
                false,
                105,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm - TSV_BBVIA30_BM1_3_BM2_3step1, 2.125*mm + TSV_BBVIA30_BM1_3_BM2_3step1, 0.5*(BM1_3High+BM2_3High)),
                logicTSV_BBVIA30_BM1_3_BM2_3,
                "TSV_BBVIA30_BM1_3_BM2_3",
                logicEnv,
                false,
                106,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm - TSV_BBVIA30_BM1_3_BM2_3step1, 2.125*mm - TSV_BBVIA30_BM1_3_BM2_3step1, 0.5*(BM1_3High+BM2_3High)),
                logicTSV_BBVIA30_BM1_3_BM2_3,
                "TSV_BBVIA30_BM1_3_BM2_3",
                logicEnv,
                false,
                107,
                true);              
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm + TSV_BBVIA30_BM1_3_BM2_3step1, -0.125*mm + TSV_BBVIA30_BM1_3_BM2_3step1, 0.5*(BM1_3High+BM2_3High)),
                logicTSV_BBVIA30_BM1_3_BM2_3,      
                "TSV_BBVIA30_BM1_3_BM2_3",
                logicEnv,
                false,
                108,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm + TSV_BBVIA30_BM1_3_BM2_3step1, -0.125*mm - TSV_BBVIA30_BM1_3_BM2_3step1, 0.5*(BM1_3High+BM2_3High)),
                logicTSV_BBVIA30_BM1_3_BM2_3,
                "TSV_BBVIA30_BM1_3_BM2_3",
                logicEnv,
                false,
                109,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm - TSV_BBVIA30_BM1_3_BM2_3step1, -0.125*mm + TSV_BBVIA30_BM1_3_BM2_3step1, 0.5*(BM1_3High+BM2_3High)),
                logicTSV_BBVIA30_BM1_3_BM2_3,
                "TSV_BBVIA30_BM1_3_BM2_3",
                logicEnv,
                false,
                110,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm - TSV_BBVIA30_BM1_3_BM2_3step1, -0.125*mm - TSV_BBVIA30_BM1_3_BM2_3step1, 0.5*(BM1_3High+BM2_3High)),
                logicTSV_BBVIA30_BM1_3_BM2_3,
                "TSV_BBVIA30_BM1_3_BM2_3",
                logicEnv,
                false,
                111,
                true);

G4double TSV_BBVIA30_BM1_3_BM2_3Xmin2 = -4.1*mm + TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Xmax2 = 0.1*mm - TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Ymin2 = -3.75*mm + TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Ymax2 = 0*mm - TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Xstep2 = (TSV_BBVIA30_BM1_3_BM2_3Xmax2 - TSV_BBVIA30_BM1_3_BM2_3Xmin2)/11;
G4double TSV_BBVIA30_BM1_3_BM2_3Ystep2 = (TSV_BBVIA30_BM1_3_BM2_3Ymax2 - TSV_BBVIA30_BM1_3_BM2_3Ymin2)/9;
for(G4int nTSV_BBVIA30_BM1_3_BM2_3X = 0; nTSV_BBVIA30_BM1_3_BM2_3X < 116; nTSV_BBVIA30_BM1_3_BM2_3X++){
    G4double x = TSV_BBVIA30_BM1_3_BM2_3Xmin2 + (nTSV_BBVIA30_BM1_3_BM2_3X%12)*TSV_BBVIA30_BM1_3_BM2_3Xstep2;
    G4double y = TSV_BBVIA30_BM1_3_BM2_3Ymin2 + (nTSV_BBVIA30_BM1_3_BM2_3X/12)*TSV_BBVIA30_BM1_3_BM2_3Ystep2;
    G4ThreeVector posTSV_BBVIA30_BM1_3_BM2_3 = G4ThreeVector(x, y, 0.5*(BM1_3High+BM2_3High));
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_3_BM2_3,
                      logicTSV_BBVIA30_BM1_3_BM2_3,
                      "TSV_BBVIA30_BM1_3_BM2_3",
                      logicEnv,
                      false,
                      nTSV_BBVIA30_BM1_3_BM2_3X + 112,
                      true);
}

G4double TSV_BBVIA30_BM1_3_BM2_3Xmin3 = -5.3*mm + TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Xmax3 = 5.3*mm - TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Ymin3 = -4.8*mm + TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Ymax3 = -3.9*mm - 5*TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Xstep3 = (TSV_BBVIA30_BM1_3_BM2_3Xmax3 - TSV_BBVIA30_BM1_3_BM2_3Xmin3)/49;
G4double TSV_BBVIA30_BM1_3_BM2_3Ystep3 = (TSV_BBVIA30_BM1_3_BM2_3Ymax3 - TSV_BBVIA30_BM1_3_BM2_3Ymin3)/2;
for(G4int nTSV_BBVIA30_BM1_3_BM2_3X = 0; nTSV_BBVIA30_BM1_3_BM2_3X < 150; nTSV_BBVIA30_BM1_3_BM2_3X++){
    G4double x = TSV_BBVIA30_BM1_3_BM2_3Xmin3 + (nTSV_BBVIA30_BM1_3_BM2_3X%50)*TSV_BBVIA30_BM1_3_BM2_3Xstep3;
    G4double y = TSV_BBVIA30_BM1_3_BM2_3Ymin3 + (nTSV_BBVIA30_BM1_3_BM2_3X/50)*TSV_BBVIA30_BM1_3_BM2_3Ystep3;
    G4ThreeVector posTSV_BBVIA30_BM1_3_BM2_3 = G4ThreeVector(x, y, 0.5*(BM1_3High+BM2_3High));
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_3_BM2_3,
                      logicTSV_BBVIA30_BM1_3_BM2_3,
                      "TSV_BBVIA30_BM1_3_BM2_3",
                      logicEnv,
                      false,
                      nTSV_BBVIA30_BM1_3_BM2_3X + 228,
                      true);
}

G4double TSV_BBVIA30_BM1_3_BM2_3Xmin4 = -5.3*mm + TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Xmax4 = 5.3*mm - TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Ymin4 = 3.9*mm + 5*TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Ymax4 = 4.8*mm - TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Xstep4 = (TSV_BBVIA30_BM1_3_BM2_3Xmax4 - TSV_BBVIA30_BM1_3_BM2_3Xmin4)/49;
G4double TSV_BBVIA30_BM1_3_BM2_3Ystep4 = (TSV_BBVIA30_BM1_3_BM2_3Ymax4 - TSV_BBVIA30_BM1_3_BM2_3Ymin4)/2;
for(G4int nTSV_BBVIA30_BM1_3_BM2_3X = 0; nTSV_BBVIA30_BM1_3_BM2_3X < 150; nTSV_BBVIA30_BM1_3_BM2_3X++){
    G4double x = TSV_BBVIA30_BM1_3_BM2_3Xmin4 + (nTSV_BBVIA30_BM1_3_BM2_3X%50)*TSV_BBVIA30_BM1_3_BM2_3Xstep4;
    G4double y = TSV_BBVIA30_BM1_3_BM2_3Ymin4 + (nTSV_BBVIA30_BM1_3_BM2_3X/50)*TSV_BBVIA30_BM1_3_BM2_3Ystep4;
    G4ThreeVector posTSV_BBVIA30_BM1_3_BM2_3 = G4ThreeVector(x, y, 0.5*(BM1_3High+BM2_3High));
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_3_BM2_3,
                      logicTSV_BBVIA30_BM1_3_BM2_3,
                      "TSV_BBVIA30_BM1_3_BM2_3",
                      logicEnv,
                      false,
                      nTSV_BBVIA30_BM1_3_BM2_3X + 378,
                      true);
}

G4double TSV_BBVIA30_BM1_3_BM2_3Xmin5 = -5.3*mm + TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Xmax5 = -4.5*mm - TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Ymin5 = -3.9*mm + 5*TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Ymax5 = 3.9*mm - 5*TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Xstep5 = (TSV_BBVIA30_BM1_3_BM2_3Xmax5 - TSV_BBVIA30_BM1_3_BM2_3Xmin5)/3;
G4double TSV_BBVIA30_BM1_3_BM2_3Ystep5 = (TSV_BBVIA30_BM1_3_BM2_3Ymax5 - TSV_BBVIA30_BM1_3_BM2_3Ymin5)/29;
for(G4int nTSV_BBVIA30_BM1_3_BM2_3X = 0; nTSV_BBVIA30_BM1_3_BM2_3X < 120; nTSV_BBVIA30_BM1_3_BM2_3X++){
    G4double x = TSV_BBVIA30_BM1_3_BM2_3Xmin5 + (nTSV_BBVIA30_BM1_3_BM2_3X%4)*TSV_BBVIA30_BM1_3_BM2_3Xstep5;
    G4double y = TSV_BBVIA30_BM1_3_BM2_3Ymin5 + (nTSV_BBVIA30_BM1_3_BM2_3X/4)*TSV_BBVIA30_BM1_3_BM2_3Ystep5;
    G4ThreeVector posTSV_BBVIA30_BM1_3_BM2_3 = G4ThreeVector(x, y, 0.5*(BM1_3High+BM2_3High));
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_3_BM2_3,
                      logicTSV_BBVIA30_BM1_3_BM2_3,
                      "TSV_BBVIA30_BM1_3_BM2_3",
                      logicEnv,
                      false,
                      nTSV_BBVIA30_BM1_3_BM2_3X + 528,
                      true);
}

G4double TSV_BBVIA30_BM1_3_BM2_3Xmin6 = 4.5*mm + TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Xmax6 = 5.3*mm - TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Ymin6 = -3.9*mm + 5*TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Ymax6 = 3.9*mm - 5*TSV_BBVIA30_BM1_3_BM2_3R;
G4double TSV_BBVIA30_BM1_3_BM2_3Xstep6 = (TSV_BBVIA30_BM1_3_BM2_3Xmax6 - TSV_BBVIA30_BM1_3_BM2_3Xmin6)/3;
G4double TSV_BBVIA30_BM1_3_BM2_3Ystep6 = (TSV_BBVIA30_BM1_3_BM2_3Ymax6 - TSV_BBVIA30_BM1_3_BM2_3Ymin6)/29;
for(G4int nTSV_BBVIA30_BM1_3_BM2_3X = 0; nTSV_BBVIA30_BM1_3_BM2_3X < 120; nTSV_BBVIA30_BM1_3_BM2_3X++){
    G4double x = TSV_BBVIA30_BM1_3_BM2_3Xmin6 + (nTSV_BBVIA30_BM1_3_BM2_3X%4)*TSV_BBVIA30_BM1_3_BM2_3Xstep6;
    G4double y = TSV_BBVIA30_BM1_3_BM2_3Ymin6 + (nTSV_BBVIA30_BM1_3_BM2_3X/4)*TSV_BBVIA30_BM1_3_BM2_3Ystep6;
    G4ThreeVector posTSV_BBVIA30_BM1_3_BM2_3 = G4ThreeVector(x, y, 0.5*(BM1_3High+BM2_3High));
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_3_BM2_3,
                      logicTSV_BBVIA30_BM1_3_BM2_3,
                      "TSV_BBVIA30_BM1_3_BM2_3",
                      logicEnv,
                      false,
                      nTSV_BBVIA30_BM1_3_BM2_3X + 648,
                      true);
}




//BPI_CORE_3,type:Dielectric,Material:FR_4                                          
G4double BPI_CORE_3High = 1.25104*mm;
G4Box* solidBPI_CORE_3 =
    new G4Box("BPI_CORE_3",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BPI_CORE_3High-BM1_3High));
G4LogicalVolume* logicBPI_CORE_3 =
    new G4LogicalVolume(solidBPI_CORE_3,
                        EpoxyResin,
                        "BPI_CORE_3");
G4ThreeVector posBPI_CORE_3 = G4ThreeVector(0,0,0.5*(BPI_CORE_3High+BM1_3High));
    new G4PVPlacement(0,
                      posBPI_CORE_3,
                      logicBPI_CORE_3,
                      "BPI_CORE_3",
                      logicEnv,
                      false,
                      0,
                      true);

//CORE2_3,type:Dielectric,Material:FR_4
G4double CORE2_3High = 1.28104*mm;  
G4Box* solidCORE2_3 =
    new G4Box("CORE2_3",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(CORE2_3High-BPI_CORE_3High));
G4LogicalVolume* logicCORE2_3 =
    new G4LogicalVolume(solidCORE2_3,
                        EpoxyResin,
                        "CORE2_3");
G4ThreeVector posCORE2_3 = G4ThreeVector(0,0,0.5*(CORE2_3High+BPI_CORE_3High));
    new G4PVPlacement(0,
                      posCORE2_3,
                      logicCORE2_3,
                      "CORE2_3",
                      logicEnv,
                      false,
                      0,
                      true);                    

//CORE1_3,type:Dielectric,Material:FR_4
G4double CORE1_3High = 1.40104*mm;
G4Box* solidCORE1_3 =
    new G4Box("CORE1_3",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(CORE1_3High-CORE2_3High));
G4LogicalVolume* logicCORE1_3 =
    new G4LogicalVolume(solidCORE1_3,
                        EpoxyResin,
                        "CORE1_3");
G4ThreeVector posCORE1_3 = G4ThreeVector(0,0,0.5*(CORE1_3High+CORE2_3High));
    new G4PVPlacement(0,
                      posCORE1_3,
                      logicCORE1_3,
                      "CORE1_3",
                      logicEnv,
                      false,
                      0,
                      true);

//FPI_CORE_3,type:Dielectric,Material:FR_4                      
G4double FPI_CORE_3High = 1.40904*mm;
G4Box* solidFPI_CORE_3 =
    new G4Box("FPI_CORE_3",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FPI_CORE_3High-CORE1_3High));
G4LogicalVolume* logicFPI_CORE_3 =
    new G4LogicalVolume(solidFPI_CORE_3,
                        EpoxyResin,
                        "FPI_CORE_3");
G4ThreeVector posFPI_CORE_3 = G4ThreeVector(0,0,0.5*(FPI_CORE_3High+CORE1_3High));
    new G4PVPlacement(0,
                      posFPI_CORE_3,
                      logicFPI_CORE_3,
                      "FPI_CORE_3",
                      logicEnv,
                      false,
                      0,
                      true);

//FM1_3,type:Metal,Material:Filled                      
G4double FM1_3High = 1.41304*mm;
G4Box* solidFM1_3 =
    new G4Box("FM1_3",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FM1_3High-FPI_CORE_3High));
G4LogicalVolume* logicFM1_3 =
    new G4LogicalVolume(solidFM1_3,
                        Filled,
                        "FM1_3");
G4ThreeVector posFM1_3 = G4ThreeVector(0,0,0.5*(FM1_3High+FPI_CORE_3High));
    new G4PVPlacement(0,
                      posFM1_3,
                      logicFM1_3,
                      "FM1_3",
                      logicEnv,
                      false,
                      0,
                      true);

//焊盘：TSV-BBVIA-30-FM1_3_BM1_3,0.01mm
G4double TSV_BBVIA_30_FM1_3_BM1_3R = 0.01*mm;
G4Tubs* solidTSV_BBVIA_30_FM1_3_BM1_3 =
    new G4Tubs("TSV_BBVIA_30_FM1_3_BM1_3",
               0,
               TSV_BBVIA_30_FM1_3_BM1_3R,
               0.5*(FM1_3High-BM1_3High),
               0.*deg,
               360.*deg);
G4LogicalVolume* logicTSV_BBVIA_30_FM1_3_BM1_3 =
    new G4LogicalVolume(solidTSV_BBVIA_30_FM1_3_BM1_3,
                        copper,
                        "TSV_BBVIA_30_FM1_3_BM1_3");
G4double TSV_BBVIA_30_FM1_3_BM1_3Xmin1 = -5.25*mm + TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Xmax1 = 5.25*mm - TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Ymin1 = -4.8*mm + TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Ymax1 = -4.05*mm - 5*TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Xstep1 = (TSV_BBVIA_30_FM1_3_BM1_3Xmax1 - TSV_BBVIA_30_FM1_3_BM1_3Xmin1)/44;
G4double TSV_BBVIA_30_FM1_3_BM1_3Ystep1 = (TSV_BBVIA_30_FM1_3_BM1_3Ymax1 - TSV_BBVIA_30_FM1_3_BM1_3Ymin1)/3;
for (int nTSV_BBVIA_30_FM1_3_BM1_3=0; nTSV_BBVIA_30_FM1_3_BM1_3<180; nTSV_BBVIA_30_FM1_3_BM1_3++){
    G4double z = 0.5 * (FM1_3High + BM1_3High);
    G4double x = TSV_BBVIA_30_FM1_3_BM1_3Xmin1 + TSV_BBVIA_30_FM1_3_BM1_3Xstep1*(nTSV_BBVIA_30_FM1_3_BM1_3%45);
    G4double y = TSV_BBVIA_30_FM1_3_BM1_3Ymin1 + TSV_BBVIA_30_FM1_3_BM1_3Ystep1*(nTSV_BBVIA_30_FM1_3_BM1_3/45);
    G4ThreeVector posTSV_BBVIA_30_FM1_3_BM1_3 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA_30_FM1_3_BM1_3,
                      logicTSV_BBVIA_30_FM1_3_BM1_3,
                      "TSV_BBVIA_30_FM1_3_BM1_3",
                      logicEnv,
                      false,
                      nTSV_BBVIA_30_FM1_3_BM1_3,
                      true);
}
G4double TSV_BBVIA_30_FM1_3_BM1_3Xmin2 = -5.25*mm + TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Xmax2 = 5.25*mm - TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Ymin2 = 4.05*mm + 5*TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Ymax2 = 4.8*mm - TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Xstep2 = (TSV_BBVIA_30_FM1_3_BM1_3Xmax2 - TSV_BBVIA_30_FM1_3_BM1_3Xmin2)/44;
G4double TSV_BBVIA_30_FM1_3_BM1_3Ystep2 = (TSV_BBVIA_30_FM1_3_BM1_3Ymax2 - TSV_BBVIA_30_FM1_3_BM1_3Ymin2)/3;
for (int nTSV_BBVIA_30_FM1_3_BM1_3=0; nTSV_BBVIA_30_FM1_3_BM1_3<180; nTSV_BBVIA_30_FM1_3_BM1_3++){
    G4double z = 0.5 * (FM1_3High + BM1_3High);
    G4double x = TSV_BBVIA_30_FM1_3_BM1_3Xmin2 + TSV_BBVIA_30_FM1_3_BM1_3Xstep2*(nTSV_BBVIA_30_FM1_3_BM1_3%45);
    G4double y = TSV_BBVIA_30_FM1_3_BM1_3Ymin2 + TSV_BBVIA_30_FM1_3_BM1_3Ystep2*(nTSV_BBVIA_30_FM1_3_BM1_3/45);
    G4ThreeVector posTSV_BBVIA_30_FM1_3_BM1_3 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA_30_FM1_3_BM1_3,
                      logicTSV_BBVIA_30_FM1_3_BM1_3,
                      "TSV_BBVIA_30_FM1_3_BM1_3",
                      logicEnv,
                      false,
                      nTSV_BBVIA_30_FM1_3_BM1_3+180,
                      true);
}
G4double TSV_BBVIA_30_FM1_3_BM1_3Xmin3 = -5.25*mm + TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Xmax3 = -4.5*mm - TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Ymin3 = -4.05*mm + 5*TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Ymax3 = 4.05*mm - 5*TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Xstep3 = (TSV_BBVIA_30_FM1_3_BM1_3Xmax3 - TSV_BBVIA_30_FM1_3_BM1_3Xmin3)/3;
G4double TSV_BBVIA_30_FM1_3_BM1_3Ystep3 = (TSV_BBVIA_30_FM1_3_BM1_3Ymax3 - TSV_BBVIA_30_FM1_3_BM1_3Ymin3)/35;
for (int nTSV_BBVIA_30_FM1_3_BM1_3=0; nTSV_BBVIA_30_FM1_3_BM1_3<143; nTSV_BBVIA_30_FM1_3_BM1_3++){
    G4double z = 0.5 * (FM1_3High + BM1_3High);
    G4double x = TSV_BBVIA_30_FM1_3_BM1_3Xmin3 + TSV_BBVIA_30_FM1_3_BM1_3Xstep3*(nTSV_BBVIA_30_FM1_3_BM1_3%4);
    G4double y = TSV_BBVIA_30_FM1_3_BM1_3Ymin3 + TSV_BBVIA_30_FM1_3_BM1_3Ystep3*(nTSV_BBVIA_30_FM1_3_BM1_3/4);
    G4ThreeVector posTSV_BBVIA_30_FM1_3_BM1_3 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA_30_FM1_3_BM1_3,
                      logicTSV_BBVIA_30_FM1_3_BM1_3,
                      "TSV_BBVIA_30_FM1_3_BM1_3",
                      logicEnv,
                      false,
                      nTSV_BBVIA_30_FM1_3_BM1_3+360,
                      true);
}
G4double TSV_BBVIA_30_FM1_3_BM1_3Xmin4 = 4.5*mm + TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Xmax4 = 5.25*mm - TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Ymin4 = -4.05*mm + 5*TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Ymax4 = 4.05*mm - 5*TSV_BBVIA_30_FM1_3_BM1_3R;
G4double TSV_BBVIA_30_FM1_3_BM1_3Xstep4 = (TSV_BBVIA_30_FM1_3_BM1_3Xmax4 - TSV_BBVIA_30_FM1_3_BM1_3Xmin4)/3;
G4double TSV_BBVIA_30_FM1_3_BM1_3Ystep4 = (TSV_BBVIA_30_FM1_3_BM1_3Ymax4 - TSV_BBVIA_30_FM1_3_BM1_3Ymin4)/35;
for (int nTSV_BBVIA_30_FM1_3_BM1_3=0; nTSV_BBVIA_30_FM1_3_BM1_3<142; nTSV_BBVIA_30_FM1_3_BM1_3++){
    G4double z = 0.5 * (FM1_3High + BM1_3High);
    G4double x = TSV_BBVIA_30_FM1_3_BM1_3Xmin4 + TSV_BBVIA_30_FM1_3_BM1_3Xstep4*(nTSV_BBVIA_30_FM1_3_BM1_3%4);
    G4double y = TSV_BBVIA_30_FM1_3_BM1_3Ymin4 + TSV_BBVIA_30_FM1_3_BM1_3Ystep4*(nTSV_BBVIA_30_FM1_3_BM1_3/4);
    G4ThreeVector posTSV_BBVIA_30_FM1_3_BM1_3 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA_30_FM1_3_BM1_3,
                      logicTSV_BBVIA_30_FM1_3_BM1_3,
                      "TSV_BBVIA_30_FM1_3_BM1_3",
                      logicEnv,
                      false,
                      nTSV_BBVIA_30_FM1_3_BM1_3+503,
                      true);
}
/////


//FP1_3,type:Dielectric,Material:FR_4
G4double FP1_3High = 1.42104*mm;  
G4Box* solidFP1_3 =
    new G4Box("FP1_3",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FP1_3High-FM1_3High));
G4LogicalVolume* logicFP1_3 =
    new G4LogicalVolume(solidFP1_3,
                        EpoxyResin,
                        "FP1_3");
G4ThreeVector posFP1_3 = G4ThreeVector(0,0,0.5*(FP1_3High+FM1_3High));
    new G4PVPlacement(0,
                      posFP1_3,
                      logicFP1_3,
                      "FP1_3",
                      logicEnv,
                      false,
                      0,
                      true);

//UBM_3,type:Metal,Material:Filled
G4double UBM_3High = 1.42504*mm;        
G4Box* solidUBM_3 =
    new G4Box("UBM_3",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UBM_3High-FP1_3High));
G4LogicalVolume* logicUBM_3 =
    new G4LogicalVolume(solidUBM_3,
                        Filled,
                        "UBM_3");
G4ThreeVector posUBM_3 = G4ThreeVector(0,0,0.5*(UBM_3High+FP1_3High));
    new G4PVPlacement(0,
                      posUBM_3,
                      logicUBM_3,
                      "UBM_3",
                      logicEnv,
                      false,
                      0,
                      true);

//UNNAMED_028,type:Dielectric,Material:Air
G4double UNNAMED_028High = 1.50504*mm;   
G4Box* solidUNNAMED_028 =
    new G4Box("UNNAMED_028",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UNNAMED_028High-UBM_3High));
G4LogicalVolume* logicUNNAMED_028 =
    new G4LogicalVolume(solidUNNAMED_028,
                        air,
                        "UNNAMED_028");
G4ThreeVector posUNNAMED_028 = G4ThreeVector(0,0,0.5*(UNNAMED_028High+UBM_3High));
    new G4PVPlacement(0,
                      posUNNAMED_028,
                      logicUNNAMED_028,
                      "UNNAMED_028",
                      logicEnv,
                      false,
                      0,
                      true);

//UBM_B_2,type:Metal,Material:Filled                     
G4double UBM_B_2High = 1.50904*mm;
G4Box* solidUBM_B_2 =
    new G4Box("UBM_B_2",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UBM_B_2High-UNNAMED_028High));
G4LogicalVolume* logicUBM_B_2 =
    new G4LogicalVolume(solidUBM_B_2,
                        Filled,
                        "UBM_B_2");
G4ThreeVector posUBM_B_2 = G4ThreeVector(0,0,0.5*(UBM_B_2High+UNNAMED_028High));
    new G4PVPlacement(0,
                      posUBM_B_2,
                      logicUBM_B_2,
                      "UBM_B_2",
                      logicEnv,
                      false,
                      0,
                      true);

//BP2_2,type:Dielectric,Material:FR_4
G4double BP2_2High = 1.71224*mm;
G4Box* solidBP2_2 =
    new G4Box("BP2_2",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BP2_2High-UBM_B_2High));
G4LogicalVolume* logicBP2_2 =
    new G4LogicalVolume(solidBP2_2,
                        EpoxyResin,
                        "BP2_2");
G4ThreeVector posBP2_2 = G4ThreeVector(0,0,0.5*(BP2_2High+UBM_B_2High));
    new G4PVPlacement(0,
                      posBP2_2,
                      logicBP2_2,
                      "BP2_2",
                      logicEnv,
                      false,
                      0,
                      true);

//BM_2_2,type:Metal,Material:Filled
G4double BM2_2High = 1.74272*mm;      
G4Box* solidBM2_2 =
    new G4Box("BM2_2",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BM2_2High-BP2_2High));
G4LogicalVolume* logicBM2_2 =
    new G4LogicalVolume(solidBM2_2,
                        Filled,
                        "BM2_2");
G4ThreeVector posBM2_2 = G4ThreeVector(0,0,0.5*(BM2_2High+BP2_2High));
    new G4PVPlacement(0,
                      posBM2_2,
                      logicBM2_2,
                      "BM2_2",
                      logicEnv,
                      false,
                      0,
                      true);
//焊盘
G4double TSVINNERBGAPAD2R = 0.03*mm;
G4Tubs* solidTSVINNERBGAPAD2 =
    new G4Tubs("TSVINNERBGAPAD2",
    0,
    TSVINNERBGAPAD2R,
    0.5*(BM2_2High-FM1_3High),
    0.*deg,
    360.*deg);
G4LogicalVolume* logicTSVINNERBGAPAD2 =
    new G4LogicalVolume(solidTSVINNERBGAPAD2,
                        copper,
                        "TSVINNERBGAPAD2");
G4int nTSVINNERBGAPAD2 = 0;                        
for(G4double x = -4*mm; x<=4*mm; x+=0.5*mm){
    for(G4double y = 0.625*mm; y<=3.625*mm; y+=0.5*mm){
        G4double z = 0.5 * (FM1_3High + BM2_2High);
        G4ThreeVector posTSVINNERBGAPAD2 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD2,
                      logicTSVINNERBGAPAD2,
                      "TSVINNERBGAPAD2",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD2,
                      true);
        nTSVINNERBGAPAD2++;
    }
}
for(G4double x = -4*mm; x<=4*mm; x+=0.5*mm){
    for(G4double y = -3.625*mm; y<=-0.125*mm; y+=0.5*mm){
        G4double z = 0.5 * (FM1_3High + BM2_2High);
        G4ThreeVector posTSVINNERBGAPAD2 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD2,
                      logicTSVINNERBGAPAD2,
                      "TSVINNERBGAPAD2",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD2,
                      true);
        nTSVINNERBGAPAD2++;
    }
}
for(G4double x = -5.125*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = 4.125*mm; y<=4.625*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_3High + BM2_2High);
        G4ThreeVector posTSVINNERBGAPAD2 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD2,
                      logicTSVINNERBGAPAD2,
                      "TSVINNERBGAPAD2",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD2,
                      true);
        nTSVINNERBGAPAD2++;
    }
}

for(G4double x = -5.125*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = -4.625*mm; y<=-4.125*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_3High + BM2_2High);
        G4ThreeVector posTSVINNERBGAPAD2 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD2,
                      logicTSVINNERBGAPAD2,
                      "TSVINNERBGAPAD2",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD2,
                      true);
        nTSVINNERBGAPAD2++;
    }
}

for(G4double x = -5.125*mm; x<=-4.625*mm; x+=0.25*mm){
    for(G4double y = -3.875*mm; y<=3.875*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_3High + BM2_2High);
        G4ThreeVector posTSVINNERBGAPAD2 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD2,
                      logicTSVINNERBGAPAD2,
                      "TSVINNERBGAPAD2",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD2,
                      true);
        nTSVINNERBGAPAD2++;
    }
}

for(G4double x = 4.625*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = -3.875*mm; y<=3.875*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_3High + BM2_2High);
        G4ThreeVector posTSVINNERBGAPAD2 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD2,
                      logicTSVINNERBGAPAD2,
                      "TSVINNERBGAPAD2",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD2,
                      true);
        nTSVINNERBGAPAD2++;
    }
}


//BP1_2,type:Dielectric,Material:FR_4
G4double BP1_2High = 1.75072*mm;
G4Box* solidBP1_2 =
    new G4Box("BP1_2",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BP1_2High-BM2_2High));
G4LogicalVolume* logicBP1_2 =
    new G4LogicalVolume(solidBP1_2,
                        EpoxyResin,
                        "BP1_2");
G4ThreeVector posBP1_2 = G4ThreeVector(0,0,0.5*(BP1_2High+BM2_2High));
    new G4PVPlacement(0,
                      posBP1_2,
                      logicBP1_2,
                      "BP1_2",
                      logicEnv,
                      false,
                      0,
                      true);

//BM1_2,type:Metal,Material:Filled
G4double BM1_2High = 1.75572*mm;
G4Box* solidBM1_2 =
    new G4Box("BM1_2",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BM1_2High-BP1_2High));
G4LogicalVolume* logicBM1_2 =
    new G4LogicalVolume(solidBM1_2,
                        Filled,
                        "BM1_2");
G4ThreeVector posBM1_2 = G4ThreeVector(0,0,0.5*(BM1_2High+BP1_2High));
    new G4PVPlacement(0,
                      posBM1_2,
                      logicBM1_2,
                      "BM1_2",
                      logicEnv,
                      false,
                      0,
                      true);

//焊盘：TSV-BBVIA30-BM1_2_BM2_2，0.015mm,812,copper
G4double TSV_BBVIA30_BM1_2_BM2_2R = 0.015*mm;    
G4Tubs* solidTSV_BBVIA30_BM1_2_BM2_2 =
    new G4Tubs("TSV_BBVIA30_BM1_2_BM2_2",
    0,
    TSV_BBVIA30_BM1_2_BM2_2R,
    0.5*(BM1_2High-BM2_2High),
    0.*deg,
    360.*deg);        
G4LogicalVolume* logicTSV_BBVIA30_BM1_2_BM2_2 =
    new G4LogicalVolume(solidTSV_BBVIA30_BM1_2_BM2_2,
                        copper,
                        "TSV_BBVIA30_BM1_2_BM2_2");
G4double TSV_BBVIA30_BM1_2_BM2_2Xmin1 = -2.1*mm + TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Xmax1 = 4.1*mm - TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Ymin1 = 0.5*mm + TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Ymax1 = 3.75*mm - TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Xstep1 = (TSV_BBVIA30_BM1_2_BM2_2Xmax1 - TSV_BBVIA30_BM1_2_BM2_2Xmin1)/18;
G4double TSV_BBVIA30_BM1_2_BM2_2Ystep1 = (TSV_BBVIA30_BM1_2_BM2_2Ymax1 - TSV_BBVIA30_BM1_2_BM2_2Ymin1)/7;
for (G4int i = 0; i < 152; i++){
    G4double x = TSV_BBVIA30_BM1_2_BM2_2Xmin1 + i%19*TSV_BBVIA30_BM1_2_BM2_2Xstep1;
    G4double y = TSV_BBVIA30_BM1_2_BM2_2Ymin1 + i/19*TSV_BBVIA30_BM1_2_BM2_2Ystep1;
    G4double z = 0.5 * (BM1_2High + BM2_2High);
    G4ThreeVector posTSV_BBVIA30_BM1_2_BM2_2 = G4ThreeVector(x, y, z);
    new G4PVPlacement(0,
                      posTSV_BBVIA30_BM1_2_BM2_2,
                      logicTSV_BBVIA30_BM1_2_BM2_2,
                      "TSV_BBVIA30_BM1_2_BM2_2",
                      logicEnv,
                      false,
                      i,
                      true);
}
G4double TSV_BBVIA30_BM1_2_BM2_2step1 = 0.06369*mm;
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -0.125*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                152,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -0.125*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                153,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -0.125*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                154,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -0.125*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                155,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -0.625*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                156,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -0.625*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                157,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -0.625*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                158,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -0.625*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                159,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -1.125*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                160,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -1.125*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                161,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -1.125*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                162,
                true);
new G4PVPlacement(0,
                G4ThreeVector(4.00*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -1.125*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                163,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -0.125*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                164,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -0.125*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                165,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -0.125*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                166,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -0.125*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                167,
                true);      
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -0.625*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                168,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -0.625*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                169,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -0.625*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                170,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -0.625*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                171,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -1.125*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                172,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -1.125*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                173,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -1.125*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                174,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -1.125*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                175,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -1.625*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                176,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -1.625*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                177,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -1.625*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                178,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-1.5*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -1.625*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                179,
                true);   
new G4PVPlacement(0,
                G4ThreeVector(-2*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -2.125*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                180,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-2*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -2.125*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                181,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(-2*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -2.125*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                182,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-2*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -2.125*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                183,
                true);     
new G4PVPlacement(0,
                G4ThreeVector(-2.5*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -2.625*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                184,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-2.5*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -2.625*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                185,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(-2.5*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -2.625*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                186,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-2.5*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -2.625*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                187,
                true);   
new G4PVPlacement(0,
                G4ThreeVector(-3.0*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -3.125*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                188,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-3.0*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -3.125*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                189,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(-3.0*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -3.125*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                190,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-3.0*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -3.125*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                191,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(-3.5*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -3.625*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                192,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-3.5*mm + TSV_BBVIA30_BM1_2_BM2_2step1, -3.625*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                193,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(-3.5*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -3.625*mm + TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                194,
                true);
new G4PVPlacement(0,
                G4ThreeVector(-3.5*mm - TSV_BBVIA30_BM1_2_BM2_2step1, -3.625*mm - TSV_BBVIA30_BM1_2_BM2_2step1, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                195,
                true);  
new G4PVPlacement(0,
                G4ThreeVector(0.125*mm, 4.81378*mm, 0.5*(BM1_2High+BM2_2High)),
                logicTSV_BBVIA30_BM1_2_BM2_2,
                "TSV_BBVIA30_BM1_2_BM2_2",
                logicEnv,
                false,
                196,
                true);  
G4double TSV_BBVIA30_BM1_2_BM2_2Xmin2 = -5.25*mm + TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Xmax2 = -4.5*mm - 5*TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Ymin2 = -4.8*mm + TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Ymax2 = 4.8*mm - TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Xstep2 = (TSV_BBVIA30_BM1_2_BM2_2Xmax2-TSV_BBVIA30_BM1_2_BM2_2Xmin2)/3;
G4double TSV_BBVIA30_BM1_2_BM2_2Ystep2 = (TSV_BBVIA30_BM1_2_BM2_2Ymax2-TSV_BBVIA30_BM1_2_BM2_2Ymin2)/40;
for (G4int i=0;i<164;i++){
    G4double x = TSV_BBVIA30_BM1_2_BM2_2Xmin2 + i%4*TSV_BBVIA30_BM1_2_BM2_2Xstep2;
    G4double y = TSV_BBVIA30_BM1_2_BM2_2Ymin2 + i/4*TSV_BBVIA30_BM1_2_BM2_2Ystep2;
    new G4PVPlacement(0,
                      G4ThreeVector(x,y,0.5*(BM1_2High+BM2_2High)),
                      logicTSV_BBVIA30_BM1_2_BM2_2,
                      "TSV_BBVIA30_BM1_2_BM2_2",
                      logicEnv,
                      false,
                      197+i,
                      true);
}
G4double TSV_BBVIA30_BM1_2_BM2_2Xmin3 = 4.5*mm + 5*TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Xmax3 = 5.25*mm - TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Ymin3 = -4.8*mm + TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Ymax3 = 4.8*mm - TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Xstep3 = (TSV_BBVIA30_BM1_2_BM2_2Xmax3-TSV_BBVIA30_BM1_2_BM2_2Xmin3)/3;
G4double TSV_BBVIA30_BM1_2_BM2_2Ystep3 = (TSV_BBVIA30_BM1_2_BM2_2Ymax3-TSV_BBVIA30_BM1_2_BM2_2Ymin3)/40;
for (G4int i=0;i<164;i++){
    G4double x = TSV_BBVIA30_BM1_2_BM2_2Xmin3 + i%4*TSV_BBVIA30_BM1_2_BM2_2Xstep3;
    G4double y = TSV_BBVIA30_BM1_2_BM2_2Ymin3 + i/4*TSV_BBVIA30_BM1_2_BM2_2Ystep3;
    new G4PVPlacement(0,
                      G4ThreeVector(x,y,0.5*(BM1_2High+BM2_2High)),
                      logicTSV_BBVIA30_BM1_2_BM2_2,
                      "TSV_BBVIA30_BM1_2_BM2_2",
                      logicEnv,
                      false,
                      361+i,
                      true);
}
G4double TSV_BBVIA30_BM1_2_BM2_2Xmin4 = -4.5*mm + 5*TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Xmax4 = 4.5*mm - 5*TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Ymin4 = -4.8*mm + TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Ymax4 = -3.9*mm - TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Xstep4 = (TSV_BBVIA30_BM1_2_BM2_2Xmax4-TSV_BBVIA30_BM1_2_BM2_2Xmin4)/35;
G4double TSV_BBVIA30_BM1_2_BM2_2Ystep4 = (TSV_BBVIA30_BM1_2_BM2_2Ymax4-TSV_BBVIA30_BM1_2_BM2_2Ymin4)/3;
for (G4int i=0;i<144;i++){
    G4double x = TSV_BBVIA30_BM1_2_BM2_2Xmin4 + i%36*TSV_BBVIA30_BM1_2_BM2_2Xstep4;
    G4double y = TSV_BBVIA30_BM1_2_BM2_2Ymin4 + i/36*TSV_BBVIA30_BM1_2_BM2_2Ystep4;
    new G4PVPlacement(0,
                      G4ThreeVector(x,y,0.5*(BM1_2High+BM2_2High)),
                      logicTSV_BBVIA30_BM1_2_BM2_2,
                      "TSV_BBVIA30_BM1_2_BM2_2",
                      logicEnv,
                      false,
                      i + 525,
                      true);
}
G4double TSV_BBVIA30_BM1_2_BM2_2Xmin5 = -4.5*mm + 5*TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Xmax5 = 4.5*mm - 5*TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Ymin5 = 3.9*mm + TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Ymax5 = 4.8*mm - TSV_BBVIA30_BM1_2_BM2_2R;
G4double TSV_BBVIA30_BM1_2_BM2_2Xstep5 = (TSV_BBVIA30_BM1_2_BM2_2Xmax5-TSV_BBVIA30_BM1_2_BM2_2Xmin5)/35;
G4double TSV_BBVIA30_BM1_2_BM2_2Ystep5 = (TSV_BBVIA30_BM1_2_BM2_2Ymax5-TSV_BBVIA30_BM1_2_BM2_2Ymin5)/3;
for (G4int i=0;i<143;i++){
    G4double x = TSV_BBVIA30_BM1_2_BM2_2Xmin5 + i%36*TSV_BBVIA30_BM1_2_BM2_2Xstep5;
    G4double y = TSV_BBVIA30_BM1_2_BM2_2Ymin5 + i/36*TSV_BBVIA30_BM1_2_BM2_2Ystep5;
    new G4PVPlacement(0,
                      G4ThreeVector(x,y,0.5*(BM1_2High+BM2_2High)),
                      logicTSV_BBVIA30_BM1_2_BM2_2,
                      "TSV_BBVIA30_BM1_2_BM2_2",
                      logicEnv,
                      false,
                      i+ 669,
                      true);
}

//BPI_CORE_2,type:Dielectric,Material:FR_4
G4double BPI_CORE_2High = 1.76372*mm;
G4Box* solidBPI_CORE_2 =
    new G4Box("BPI_CORE_2",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BPI_CORE_2High-BM1_2High));
G4LogicalVolume* logicBPI_CORE_2 =
    new G4LogicalVolume(solidBPI_CORE_2,
                        EpoxyResin,
                        "BPI_CORE_2");
G4ThreeVector posBPI_CORE_2 = G4ThreeVector(0,0,0.5*(BPI_CORE_2High+BM1_2High));
    new G4PVPlacement(0,
                      posBPI_CORE_2,
                      logicBPI_CORE_2,
                      "BPI_CORE_2",
                      logicEnv,
                      false,
                      0,
                      true);

//CORE2_2,type:Dielectric,Material:FR_4
G4double CORE2_2High = 1.79372*mm;
G4Box* solidCORE2_2 =
    new G4Box("CORE2_2",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(CORE2_2High-BPI_CORE_2High));
G4LogicalVolume* logicCORE2_2 =
    new G4LogicalVolume(solidCORE2_2,
                        EpoxyResin,
                        "CORE2_2");
G4ThreeVector posCORE2_2 = G4ThreeVector(0,0,0.5*(CORE2_2High+BPI_CORE_2High));
    new G4PVPlacement(0,
                      posCORE2_2,
                      logicCORE2_2,
                      "CORE2_2",
                      logicEnv,
                      false,
                      0,
                      true);

//CORE1_2,type:Dielectric,Material:FR_4
G4double CORE1_2High = 1.91372*mm;
G4Box* solidCORE1_2 =
    new G4Box("CORE1_2",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(CORE1_2High-CORE2_2High));
G4LogicalVolume* logicCORE1_2 =
    new G4LogicalVolume(solidCORE1_2,
                        EpoxyResin,
                        "CORE1_2");
G4ThreeVector posCORE1_2 = G4ThreeVector(0,0,0.5*(CORE1_2High+CORE2_2High));
    new G4PVPlacement(0,
                      posCORE1_2,
                      logicCORE1_2,
                      "CORE1_2",
                      logicEnv,
                      false,
                      0,
                      true);

//FPI_CORE_2,type:Dielectric,Material:FR_4
G4double FPI_CORE_2High = 1.92172*mm; 
G4Box* solidFPI_CORE_2 =
    new G4Box("FPI_CORE_2",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FPI_CORE_2High-CORE1_2High));
G4LogicalVolume* logicFPI_CORE_2 =
    new G4LogicalVolume(solidFPI_CORE_2,
                        EpoxyResin,
                        "FPI_CORE_2");
G4ThreeVector posFPI_CORE_2 = G4ThreeVector(0,0,0.5*(FPI_CORE_2High+CORE1_2High));
    new G4PVPlacement(0,
                      posFPI_CORE_2,
                      logicFPI_CORE_2,
                      "FPI_CORE_2",
                      logicEnv,
                      false,
                      0,
                      true);

//FM1_2,type:Metal,Material:Filled
G4double FM1_2High = 1.92972*mm;
G4Box* solidFM1_2 =
    new G4Box("FM1_2",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FM1_2High-FPI_CORE_2High));
G4LogicalVolume* logicFM1_2 =
    new G4LogicalVolume(solidFM1_2,
                        Filled,
                        "FM1_2");
G4ThreeVector posFM1_2 = G4ThreeVector(0,0,0.5*(FM1_2High+FPI_CORE_2High));
    new G4PVPlacement(0,
                      posFM1_2,
                      logicFM1_2,
                      "FM1_2",
                      logicEnv,
                      false,
                      0,
                      true);

//TSV-BBVIA-30-FM1_2_BM1_2，0.01mm,copper,焊盘
G4double TSV_BBVIA30_FM1_2_BM1_2R = 0.01*mm;
G4Tubs* solidTSV_BBVIA30_FM1_2_BM1_2 =
      new G4Tubs("TSV_BBVIA30_FM1_2_BM1_2",
                0,
                TSV_BBVIA30_FM1_2_BM1_2R,
                0.5*(FM1_2High-BM1_2High),
                0.*deg,
                360.*deg);                       
G4LogicalVolume* logicTSV_BBVIA30_FM1_2_BM1_2 =
    new G4LogicalVolume(solidTSV_BBVIA30_FM1_2_BM1_2,
                        copper,
                        "TSV_BBVIA30_FM1_2_BM1_2");
G4double TSV_BBVIA30_FM1_2_BM1_2Xmin1 = -5.25*mm + TSV_BBVIA30_FM1_2_BM1_2R;                        
G4double TSV_BBVIA30_FM1_2_BM1_2Xmax1 = -4.5*mm - 5*TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Ymin1 = -4.8*mm + TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Ymax1 = 4.8*mm - TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Xstep1 = (TSV_BBVIA30_FM1_2_BM1_2Xmax1 - TSV_BBVIA30_FM1_2_BM1_2Xmin1)/3;
G4double TSV_BBVIA30_FM1_2_BM1_2Ystep1 = (TSV_BBVIA30_FM1_2_BM1_2Ymax1 - TSV_BBVIA30_FM1_2_BM1_2Ymin1)/39;
for (G4int i = 0; i < 160;i ++){
    G4double x = TSV_BBVIA30_FM1_2_BM1_2Xmin1 + i%4*TSV_BBVIA30_FM1_2_BM1_2Xstep1;
    G4double y = TSV_BBVIA30_FM1_2_BM1_2Ymin1 + i/4*TSV_BBVIA30_FM1_2_BM1_2Ystep1;
    new G4PVPlacement(0,
                      G4ThreeVector(x,y,0.5*(FM1_2High+BM1_2High)),
                      logicTSV_BBVIA30_FM1_2_BM1_2,
                      "TSV_BBVIA30_FM1_2_BM1_2",
                      logicEnv,
                      false,
                      i,
                      true);
}
G4double TSV_BBVIA30_FM1_2_BM1_2Xmin2 = 4.5*mm - 5*TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Xmax2 = 5.25*mm - TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Ymin2 = -4.8*mm + TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Ymax2 = 4.8*mm - TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Xstep2 = (TSV_BBVIA30_FM1_2_BM1_2Xmax2 - TSV_BBVIA30_FM1_2_BM1_2Xmin2)/3;
G4double TSV_BBVIA30_FM1_2_BM1_2Ystep2 = (TSV_BBVIA30_FM1_2_BM1_2Ymax2 - TSV_BBVIA30_FM1_2_BM1_2Ymin2)/39;
for (G4int i = 0; i < 160;i ++){
    G4double x = TSV_BBVIA30_FM1_2_BM1_2Xmin2 + i%4*TSV_BBVIA30_FM1_2_BM1_2Xstep2;
    G4double y = TSV_BBVIA30_FM1_2_BM1_2Ymin2 + i/4*TSV_BBVIA30_FM1_2_BM1_2Ystep2;
    new G4PVPlacement(0,
                      G4ThreeVector(x,y,0.5*(FM1_2High+BM1_2High)),
                      logicTSV_BBVIA30_FM1_2_BM1_2,
                      "TSV_BBVIA30_FM1_2_BM1_2",
                      logicEnv,
                      false,
                      i+160,
                      true);
}
G4double TSV_BBVIA30_FM1_2_BM1_2Xmin3 = -4.5*mm + 5*TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Xmax3 = 4.5*mm + 5*TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Ymin3 = -4.8*mm + TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Ymax3 = -4.05*mm - TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Xstep3 = (TSV_BBVIA30_FM1_2_BM1_2Xmax3 - TSV_BBVIA30_FM1_2_BM1_2Xmin3)/39;
G4double TSV_BBVIA30_FM1_2_BM1_2Ystep3 = (TSV_BBVIA30_FM1_2_BM1_2Ymax3 - TSV_BBVIA30_FM1_2_BM1_2Ymin3)/3;
for (G4int i = 0; i < 160;i ++){
    G4double x = TSV_BBVIA30_FM1_2_BM1_2Xmin3 + i%40*TSV_BBVIA30_FM1_2_BM1_2Xstep3;
    G4double y = TSV_BBVIA30_FM1_2_BM1_2Ymin3 + i/40*TSV_BBVIA30_FM1_2_BM1_2Ystep3;
    new G4PVPlacement(0,
                      G4ThreeVector(x,y,0.5*(FM1_2High+BM1_2High)),
                      logicTSV_BBVIA30_FM1_2_BM1_2,
                      "TSV_BBVIA30_FM1_2_BM1_2",
                      logicEnv,
                      false,
                      i+320,
                      true);
}
G4double TSV_BBVIA30_FM1_2_BM1_2Xmin4 = -4.5*mm + 5*TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Xmax4 = 4.5*mm + 5*TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Ymin4 = 4.05*mm + TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Ymax4 = 4.8*mm - TSV_BBVIA30_FM1_2_BM1_2R;
G4double TSV_BBVIA30_FM1_2_BM1_2Xstep4 = (TSV_BBVIA30_FM1_2_BM1_2Xmax4 - TSV_BBVIA30_FM1_2_BM1_2Xmin4)/39;
G4double TSV_BBVIA30_FM1_2_BM1_2Ystep4 = (TSV_BBVIA30_FM1_2_BM1_2Ymax4 - TSV_BBVIA30_FM1_2_BM1_2Ymin4)/3;
for (G4int i = 0; i < 159;i ++){
    G4double x = TSV_BBVIA30_FM1_2_BM1_2Xmin4 + i%40*TSV_BBVIA30_FM1_2_BM1_2Xstep4;
    G4double y = TSV_BBVIA30_FM1_2_BM1_2Ymin4 + i/40*TSV_BBVIA30_FM1_2_BM1_2Ystep4;
    new G4PVPlacement(0,
                      G4ThreeVector(x,y,0.5*(FM1_2High+BM1_2High)),
                      logicTSV_BBVIA30_FM1_2_BM1_2,
                      "TSV_BBVIA30_FM1_2_BM1_2",
                      logicEnv,
                      false,
                      i+480,
                      true);
}


//FP1_2,type:Dielectric,Material:FR_4
G4double FP1_2High = 1.93372*mm;      
G4Box* solidFP1_2 =
    new G4Box("FP1_2",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FP1_2High-FM1_2High));
G4LogicalVolume* logicFP1_2 =
    new G4LogicalVolume(solidFP1_2,
                        EpoxyResin,
                        "FP1_2");
G4ThreeVector posFP1_2 = G4ThreeVector(0,0,0.5*(FP1_2High+FM1_2High));
    new G4PVPlacement(0,
                      posFP1_2,
                      logicFP1_2,
                      "FP1_2",
                      logicEnv,
                      false,
                      0,
                      true);  

//UBM_2,type:Metal,Material:Filled
G4double UBM_2High = 1.93772*mm;
G4Box* solidUBM_2 =
    new G4Box("UBM_2",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UBM_2High-FP1_2High));
G4LogicalVolume* logicUBM_2 =
    new G4LogicalVolume(solidUBM_2,
                        Filled,
                        "UBM_2");
G4ThreeVector posUBM_2 = G4ThreeVector(0,0,0.5*(UBM_2High+FP1_2High));
    new G4PVPlacement(0,
                      posUBM_2,
                      logicUBM_2,
                      "UBM_2",
                      logicEnv,
                      false,
                      0,
                      true);

//UNNAMED_015,type:Dielectric,Material:Air
G4double UNNAMED_015High = 2.01772*mm;
G4Box* solidUNNAMED_015 =
    new G4Box("UNNAMED_015",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UNNAMED_015High-UBM_2High));
G4LogicalVolume* logicUNNAMED_015 =
    new G4LogicalVolume(solidUNNAMED_015,
                        air,
                        "UNNAMED_015");
G4ThreeVector posUNNAMED_015 = G4ThreeVector(0,0,0.5*(UNNAMED_015High+UBM_2High));
    new G4PVPlacement(0,
                      posUNNAMED_015,
                      logicUNNAMED_015,
                      "UNNAMED_015",
                      logicEnv,
                      false,
                      0,
                      true);

//UBM_B_1,type:Metal,Material:Filled
G4double UBM_B_1High = 2.02172*mm;
G4Box* solidUBM_B_1 =
    new G4Box("UBM_B_1",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UBM_B_1High-UNNAMED_015High));
G4LogicalVolume* logicUBM_B_1 =
    new G4LogicalVolume(solidUBM_B_1,
                        Filled,
                        "UBM_B_1");
G4ThreeVector posUBM_B_1 = G4ThreeVector(0,0,0.5*(UBM_B_1High+UNNAMED_015High));
    new G4PVPlacement(0,
                      posUBM_B_1,
                      logicUBM_B_1,
                      "UBM_B_1",
                      logicEnv,
                      false,
                      0,
                      true);
                    

//BP2_1,type:Dielectric,Material:FR_4
G4double BP2_1High = 2.22492*mm; 
G4Box* solidBP2_1 =
    new G4Box("BP2_1",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BP2_1High-UBM_B_1High));
G4LogicalVolume* logicBP2_1 =
    new G4LogicalVolume(solidBP2_1,
                        EpoxyResin,
                        "BP2_1");
G4ThreeVector posBP2_1 = G4ThreeVector(0,0,0.5*(BP2_1High+UBM_B_1High));
    new G4PVPlacement(0,
                      posBP2_1,
                      logicBP2_1,
                      "BP2_1",
                      logicEnv,
                      false,
                      0,
                      true);

//BM2_1,type:Metal,Material:Filled
G4double BM2_1High = 2.2554*mm;     
G4Box* solidBM2_1 =
    new G4Box("BM2_1",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BM2_1High-BP2_1High));
G4LogicalVolume* logicBM2_1 =
    new G4LogicalVolume(solidBM2_1,
                        Filled,
                        "BM2_1");
G4ThreeVector posBM2_1 = G4ThreeVector(0,0,0.5*(BM2_1High+BP2_1High));
    new G4PVPlacement(0,
                      posBM2_1,
                      logicBM2_1,
                      "BM2_1",
                      logicEnv,
                      false,
                      0,
                      true);
//焊盘
G4double TSVINNERBGAPAD1R = 0.03*mm;
G4Tubs* solidTSVINNERBGAPAD1 =
    new G4Tubs("TSVINNERBGAPAD1",
    0,
    TSVINNERBGAPAD1R,
    0.5*(BM2_1High-FM1_2High),
    0.*deg,
    360.*deg);
G4LogicalVolume* logicTSVINNERBGAPAD1 =
    new G4LogicalVolume(solidTSVINNERBGAPAD1,
                        copper,
                        "TSVINNERBGAPAD1");
G4int nTSVINNERBGAPAD1 = 0;                        
for(G4double x = -4*mm; x<=4*mm; x+=0.5*mm){
    for(G4double y = 0.625*mm; y<=3.625*mm; y+=0.5*mm){
        G4double z = 0.5 * (FM1_2High + BM2_1High);
        G4ThreeVector posTSVINNERBGAPAD1 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD1,
                      logicTSVINNERBGAPAD1,
                      "TSVINNERBGAPAD1",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD1,
                      true);
        nTSVINNERBGAPAD1++;
    }
}
for(G4double x = -4*mm; x<=4*mm; x+=0.5*mm){
    for(G4double y = -3.625*mm; y<=-0.125*mm; y+=0.5*mm){
        G4double z = 0.5 * (FM1_2High + BM2_1High);
        G4ThreeVector posTSVINNERBGAPAD1 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD1,
                      logicTSVINNERBGAPAD1,
                      "TSVINNERBGAPAD1",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD1,
                      true);
        nTSVINNERBGAPAD1++;
    }
}
for(G4double x = -5.125*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = 4.125*mm; y<=4.625*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_2High + BM2_1High);
        G4ThreeVector posTSVINNERBGAPAD1 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD1,
                      logicTSVINNERBGAPAD1,
                      "TSVINNERBGAPAD1",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD1,
                      true);
        nTSVINNERBGAPAD1++;
    }
}

for(G4double x = -5.125*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = -4.625*mm; y<=-4.125*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_2High + BM2_1High);
        G4ThreeVector posTSVINNERBGAPAD1 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD1,
                      logicTSVINNERBGAPAD1,
                      "TSVINNERBGAPAD1",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD1,
                      true);
        nTSVINNERBGAPAD1++;
    }
}

for(G4double x = -5.125*mm; x<=-4.625*mm; x+=0.25*mm){
    for(G4double y = -3.875*mm; y<=3.875*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_2High + BM2_1High);
        G4ThreeVector posTSVINNERBGAPAD1 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD1,
                      logicTSVINNERBGAPAD1,
                      "TSVINNERBGAPAD1",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD1,
                      true);
        nTSVINNERBGAPAD1++;
    }
}

for(G4double x = 4.625*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = -3.875*mm; y<=3.875*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_2High + BM2_1High);
        G4ThreeVector posTSVINNERBGAPAD1 = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVINNERBGAPAD1,
                      logicTSVINNERBGAPAD1,
                      "TSVINNERBGAPAD1",
                      logicEnv,
                      false,
                      nTSVINNERBGAPAD1,
                      true);
        nTSVINNERBGAPAD1++;
    }
}                      

//BP1_1,type:Dielectric,Material:Polyimide
G4double BP1_1High = 2.2634*mm;
G4Box* solidBP1_1 =
    new G4Box("BP1_1",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BP1_1High-BM2_1High));
G4LogicalVolume* logicBP1_1 =
    new G4LogicalVolume(solidBP1_1,
                        PolyimideResin,
                        "BP1_1");
G4ThreeVector posBP1_1 = G4ThreeVector(0,0,0.5*(BP1_1High+BM2_1High));
    new G4PVPlacement(0,
                      posBP1_1,
                      logicBP1_1,
                      "BP1_1",
                      logicEnv,
                      false,
                      0,
                      true);

//BM1_1,type:Metal,Material:Filled
//G4double BP1_1High = 2.2634*mm;
G4double BM1_1High = 2.2674*mm;   
G4Box* solidBM1_1 =
    new G4Box("BM1_1",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BM1_1High-BP1_1High));
G4LogicalVolume* logicBM1_1 =
    new G4LogicalVolume(solidBM1_1,
                        Filled,
                        "BM1_1");
G4ThreeVector posBM1_1 = G4ThreeVector(0,0,0.5*(BM1_1High+BP1_1High));
    new G4PVPlacement(0,
                      posBM1_1,
                      logicBM1_1,
                      "BM1_1",
                      logicEnv,
                      false,
                      0,
                      true); 

//焊盘
G4double TSV_BBVIA30_BM1_1_BM2_1R = 0.015*mm;
G4Tubs* solidTSV_BBVIA30_BM1_1_BM2_1 =
    new G4Tubs("TSV_BBVIA30_BM1_1_BM2_1",
    0,
    TSV_BBVIA30_BM1_1_BM2_1R,
    0.5*(BM1_1High-BM2_1High),
    0.*deg,
    360.*deg);
G4LogicalVolume* logicTSV_BBVIA30_BM1_1_BM2_1 =
    new G4LogicalVolume(solidTSV_BBVIA30_BM1_1_BM2_1,
                        copper,
                        "TSV_BBVIA30_BM1_1_BM2_1");
G4double TSV_BBVIA30_BM1_1_BM2_1step1 = 0.0638*mm;
G4int nTSV_BBVIA30_BM1_1_BM2_1 = 0;
for(G4double x = -3.5*mm, y = 3.625*mm; x<=-1.5*mm, y>=1.625*mm;x += 0.5*mm ,y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}          

for(G4double x = -1.0*mm, y = 3.625*mm; x<=1.0*mm, y>=1.625*mm;x += 0.5*mm ,y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}        

for(G4double x = -1.5*mm, y = 1.125*mm;y>=0.625*mm;y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
for(G4double x = 3.5*mm, y = 3.625*mm;x<=4.0*mm, y>=3.125*mm;x+=0.5*mm, y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}

for(G4double x = 3.0*mm, y = -0.125*mm; y>=-3.625*mm;y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
for(G4double x = 1.5*mm, y = -0.125*mm; y>=-3.625*mm; y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
for(G4double x = 1.0*mm, y = -0.125*mm; y>=-3.125*mm; y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
for(G4double x = 0.5*mm, y = -0.125*mm; y>=-2.625*mm; y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
for(G4double x = 0.0*mm, y = -0.125*mm; y>=-2.125*mm; y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
for(G4double x = -2.0*mm, y = -0.125*mm; y>=-1.125*mm; y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
for(G4double x = 0.5*mm, y = -3.625*mm; y>=-3.625*mm; y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
for(G4double x = 0.0*mm, y = -3.125*mm; x>=-0.5*mm, y>=-3.625*mm; x -= 0.5*mm, y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
for(G4double x = -0.5*mm, y = -2.625*mm; x>=-1.5*mm, y>=-3.625*mm; x -= 0.5*mm, y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
for(G4double x = -0.5*mm, y = -0.125*mm; y>=-0.125*mm; y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
for(G4double x = -2.5*mm, y = -1.625*mm; x>=-4.0*mm,y>=-3.125*mm; x -= 0.5*mm,y -= 0.5*mm){
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1,y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x + TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);  
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y + TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
  new G4PVPlacement(0,
                  G4ThreeVector(x - TSV_BBVIA30_BM1_1_BM2_1step1, y - TSV_BBVIA30_BM1_1_BM2_1step1, 0.5*(BM1_1High+BM2_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
G4double TSV_BBVIA30_BM1_1_BM2_1Xmin1 = -5.3*mm + TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Xmax1 = 5.3*mm - TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Ymin1 = -4.8*mm + TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Ymax1 = -4.0*mm - 5*TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Xstep1 = (TSV_BBVIA30_BM1_1_BM2_1Xmax1-TSV_BBVIA30_BM1_1_BM2_1Xmin1)/39;
G4double TSV_BBVIA30_BM1_1_BM2_1Ystep1 = (TSV_BBVIA30_BM1_1_BM2_1Ymax1-TSV_BBVIA30_BM1_1_BM2_1Ymin1)/3;
for (G4int i = 0; i < 160; i++){
  G4double x = TSV_BBVIA30_BM1_1_BM2_1Xmin1 + (i%40)*TSV_BBVIA30_BM1_1_BM2_1Xstep1;
  G4double y = TSV_BBVIA30_BM1_1_BM2_1Ymin1 + (i/40)*TSV_BBVIA30_BM1_1_BM2_1Ystep1;
  new G4PVPlacement(0,
                  G4ThreeVector(x,y,0.5*(BM2_1High+BM1_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
G4double TSV_BBVIA30_BM1_1_BM2_1Xmin2 = -5.3*mm + TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Xmax2 = 5.3*mm - TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Ymin2 = 4.0*mm + 5*TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Ymax2 = 4.8*mm - TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Xstep2 = (TSV_BBVIA30_BM1_1_BM2_1Xmax2-TSV_BBVIA30_BM1_1_BM2_1Xmin2)/39;
G4double TSV_BBVIA30_BM1_1_BM2_1Ystep2 = (TSV_BBVIA30_BM1_1_BM2_1Ymax2-TSV_BBVIA30_BM1_1_BM2_1Ymin2)/3;
for (G4int i = 0; i < 159; i++){
  G4double x = TSV_BBVIA30_BM1_1_BM2_1Xmin2 + (i%40)*TSV_BBVIA30_BM1_1_BM2_1Xstep2;
  G4double y = TSV_BBVIA30_BM1_1_BM2_1Ymin2 + (i/40)*TSV_BBVIA30_BM1_1_BM2_1Ystep2;
  new G4PVPlacement(0,
                  G4ThreeVector(x,y,0.5*(BM2_1High+BM1_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
G4double TSV_BBVIA30_BM1_1_BM2_1Xmin3 = -5.3*mm + TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Xmax3 = -4.5*mm - TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Ymin3 = -4.0*mm + 5*TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Ymax3 = 4.0*mm - 5*TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Xstep3 = (TSV_BBVIA30_BM1_1_BM2_1Xmax3-TSV_BBVIA30_BM1_1_BM2_1Xmin3)/3;
G4double TSV_BBVIA30_BM1_1_BM2_1Ystep3 = (TSV_BBVIA30_BM1_1_BM2_1Ymax3-TSV_BBVIA30_BM1_1_BM2_1Ymin3)/29;
for (G4int i = 0; i < 120; i++){
  G4double x = TSV_BBVIA30_BM1_1_BM2_1Xmin3 + (i%4)*TSV_BBVIA30_BM1_1_BM2_1Xstep3;
  G4double y = TSV_BBVIA30_BM1_1_BM2_1Ymin3 + (i/4)*TSV_BBVIA30_BM1_1_BM2_1Ystep3;
  new G4PVPlacement(0,
                  G4ThreeVector(x,y,0.5*(BM2_1High+BM1_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
G4double TSV_BBVIA30_BM1_1_BM2_1Xmin4 = 4.5*mm + TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Xmax4 = 5.3*mm - TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Ymin4 = -4.0*mm + 5*TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Ymax4 = 4.0*mm - 5*TSV_BBVIA30_BM1_1_BM2_1R;
G4double TSV_BBVIA30_BM1_1_BM2_1Xstep4 = (TSV_BBVIA30_BM1_1_BM2_1Xmax4-TSV_BBVIA30_BM1_1_BM2_1Xmin4)/3;
G4double TSV_BBVIA30_BM1_1_BM2_1Ystep4 = (TSV_BBVIA30_BM1_1_BM2_1Ymax4-TSV_BBVIA30_BM1_1_BM2_1Ymin4)/29;
for (G4int i = 0; i < 120; i++){
  G4double x = TSV_BBVIA30_BM1_1_BM2_1Xmin4 + (i%4)*TSV_BBVIA30_BM1_1_BM2_1Xstep4;
  G4double y = TSV_BBVIA30_BM1_1_BM2_1Ymin4 + (i/4)*TSV_BBVIA30_BM1_1_BM2_1Ystep4;
  new G4PVPlacement(0,
                  G4ThreeVector(x,y,0.5*(BM2_1High+BM1_1High)),
                  logicTSV_BBVIA30_BM1_1_BM2_1,
                  "TSV_BBVIA30_BM1_1_BM2_1",
                  logicEnv,
                  false,
                  nTSV_BBVIA30_BM1_1_BM2_1++,
                  true);
}
                

//BPI_CORE_1,type:Dielectric,Material:SiliconDioxide                                        
G4double BPI_CORE_1High = 2.2754*mm;
G4Box* solidBPI_CORE_1 =
    new G4Box("BPI_CORE_1",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(BPI_CORE_1High-BM1_1High));
G4LogicalVolume* logicBPI_CORE_1 =
    new G4LogicalVolume(solidBPI_CORE_1,
                        SiliconDioxide,
                        "BPI_CORE_1");
G4ThreeVector posBPI_CORE_1 = G4ThreeVector(0,0,0.5*(BPI_CORE_1High+BM1_1High));
    new G4PVPlacement(0,
                      posBPI_CORE_1,
                      logicBPI_CORE_1,
                      "BPI_CORE_1",
                      logicEnv,
                      false,
                      0,
                      true);

//CORE2_1,type:Dielectric,Material:Silicon
G4double CORE2_1High = 2.3054*mm;
G4Box* solidCORE2_1 =
    new G4Box("CORE2_1",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(CORE2_1High-BPI_CORE_1High));
G4LogicalVolume* logicCORE2_1 =
    new G4LogicalVolume(solidCORE2_1,
                        Silicon,
                        "CORE2_1");
G4ThreeVector posCORE2_1 = G4ThreeVector(0,0,0.5*(CORE2_1High+BPI_CORE_1High));
    new G4PVPlacement(0,
                      posCORE2_1,
                      logicCORE2_1,
                      "CORE2_1",
                      logicEnv,
                      false,
                      0,
                      true);

//CORE1_1,type:Dielectric,Material:Silicon
G4double CORE1_1High = 2.4254*mm;
G4Box* solidCORE1_1 =
    new G4Box("CORE1_1",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(CORE1_1High-CORE2_1High));
G4LogicalVolume* logicCORE1_1 =
    new G4LogicalVolume(solidCORE1_1,
                        Silicon,
                        "CORE1_1");
G4ThreeVector posCORE1_1 = G4ThreeVector(0,0,0.5*(CORE1_1High+CORE2_1High));
    new G4PVPlacement(0,
                      posCORE1_1,
                      logicCORE1_1,
                      "CORE1_1",
                      logicEnv,
                      false,
                      0,
                      true);

//FPI_CORE_1,type:Dielectric,Material:SiliconDioxide
G4double FPI_CORE_1High = 2.4334*mm;
G4Box* solidFPI_CORE_1 =
    new G4Box("FPI_CORE_1",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FPI_CORE_1High-CORE1_1High));
G4LogicalVolume* logicFPI_CORE_1 =
    new G4LogicalVolume(solidFPI_CORE_1,
                        SiliconDioxide,
                        "FPI_CORE_1");
G4ThreeVector posFPI_CORE_1 = G4ThreeVector(0,0,0.5*(FPI_CORE_1High+CORE1_1High));
    new G4PVPlacement(0,
                      posFPI_CORE_1,
                      logicFPI_CORE_1,
                      "FPI_CORE_1",
                      logicEnv,
                      false,
                      0,
                      true);
//FM1_1,type:Metal,Material:Filled
G4double FM1_1High = 2.4374*mm;
G4Box* solidFM1_1 =
    new G4Box("FM1_1",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FM1_1High-FPI_CORE_1High));
G4LogicalVolume* logicFM1_1 =
    new G4LogicalVolume(solidFM1_1,
                        Filled,
                        "FM1_1");
G4ThreeVector posFM1_1 = G4ThreeVector(0,0,0.5*(FM1_1High+FPI_CORE_1High));
    new G4PVPlacement(0,
                      posFM1_1,
                      logicFM1_1,
                      "FM1_1",
                      logicEnv,
                      false,
                      0,
                      true);

//焊盘
G4double TSV_BBVIA_30_FM1_1_BM1_1R = 0.01*mm;
G4Tubs* solidTSV_BBVIA_30_FM1_1_BM1_1 =
    new G4Tubs("TSV_BBVIA_30_FM1_1_BM1_1",
               0,
               TSV_BBVIA_30_FM1_1_BM1_1R,
               0.5*(FM1_1High-BM1_1High),
               0.*deg,
               360.*deg);
G4LogicalVolume* logicTSV_BBVIA_30_FM1_1_BM1_1 =
    new G4LogicalVolume(solidTSV_BBVIA_30_FM1_1_BM1_1,
                        copper,
                        "TSV_BBVIA_30_FM1_1_BM1_1");
G4int nTSV_BBVIA_30_FM1_1_BM1_1 = 0;
G4double TSV_BBVIA_30_FM1_1_BM1_1Xmin1 = -5.25*mm + TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Xmax1 = 5.25*mm - TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Ymin1 = -4.8*mm + TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Ymax1 = -4.05*mm - 5*TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Xstep1 = (TSV_BBVIA_30_FM1_1_BM1_1Xmax1 - TSV_BBVIA_30_FM1_1_BM1_1Xmin1)/42;
G4double TSV_BBVIA_30_FM1_1_BM1_1Ystep1 = (TSV_BBVIA_30_FM1_1_BM1_1Ymax1 - TSV_BBVIA_30_FM1_1_BM1_1Ymin1)/3;
for (G4int i=0; i<171; i++){
    G4double x = TSV_BBVIA_30_FM1_1_BM1_1Xmin1 + TSV_BBVIA_30_FM1_1_BM1_1Xstep1*(i%43);
    G4double y = TSV_BBVIA_30_FM1_1_BM1_1Ymin1 + TSV_BBVIA_30_FM1_1_BM1_1Ystep1*(i/43);
    new G4PVPlacement(0,
                      G4ThreeVector(x, y, 0.5 * (FM1_1High + BM1_1High)),
                      logicTSV_BBVIA_30_FM1_1_BM1_1,
                      "TSV_BBVIA_30_FM1_1_BM1_1",
                      logicEnv,
                      false,
                      nTSV_BBVIA_30_FM1_1_BM1_1++,
                      true);
}
G4double TSV_BBVIA_30_FM1_1_BM1_1Xmin2 = -5.25*mm + TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Xmax2 = 5.25*mm - TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Ymin2 = 4.05*mm + 5*TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Ymax2 = 4.8*mm - TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Xstep2 = (TSV_BBVIA_30_FM1_1_BM1_1Xmax2 - TSV_BBVIA_30_FM1_1_BM1_1Xmin2)/42;
G4double TSV_BBVIA_30_FM1_1_BM1_1Ystep2 = (TSV_BBVIA_30_FM1_1_BM1_1Ymax2 - TSV_BBVIA_30_FM1_1_BM1_1Ymin2)/3;
for (G4int i=0; i<171; i++){
    G4double x = TSV_BBVIA_30_FM1_1_BM1_1Xmin2 + TSV_BBVIA_30_FM1_1_BM1_1Xstep2*(i%43);
    G4double y = TSV_BBVIA_30_FM1_1_BM1_1Ymin2 + TSV_BBVIA_30_FM1_1_BM1_1Ystep2*(i/43);
    new G4PVPlacement(0,
                      G4ThreeVector(x, y, 0.5 * (FM1_1High + BM1_1High)),
                      logicTSV_BBVIA_30_FM1_1_BM1_1,
                      "TSV_BBVIA_30_FM1_1_BM1_1",
                      logicEnv,
                      false,
                      nTSV_BBVIA_30_FM1_1_BM1_1++,
                      true);
}
G4double TSV_BBVIA_30_FM1_1_BM1_1Xmin3 = -5.25*mm + TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Xmax3 = -4.5*mm - TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Ymin3 = -4.05*mm + 5*TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Ymax3 = 4.05*mm - 5*TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Xstep3 = (TSV_BBVIA_30_FM1_1_BM1_1Xmax3 - TSV_BBVIA_30_FM1_1_BM1_1Xmin3)/3;
G4double TSV_BBVIA_30_FM1_1_BM1_1Ystep3 = (TSV_BBVIA_30_FM1_1_BM1_1Ymax3 - TSV_BBVIA_30_FM1_1_BM1_1Ymin3)/32;
for (G4int i=0; i<131; i++){
    G4double x = TSV_BBVIA_30_FM1_1_BM1_1Xmin3 + TSV_BBVIA_30_FM1_1_BM1_1Xstep3*(i%4);
    G4double y = TSV_BBVIA_30_FM1_1_BM1_1Ymin3 + TSV_BBVIA_30_FM1_1_BM1_1Ystep3*(i/4);
    new G4PVPlacement(0,
                      G4ThreeVector(x, y, 0.5 * (FM1_1High + BM1_1High)),
                      logicTSV_BBVIA_30_FM1_1_BM1_1,
                      "TSV_BBVIA_30_FM1_1_BM1_1",
                      logicEnv,
                      false,
                      nTSV_BBVIA_30_FM1_1_BM1_1++,
                      true);
}
G4double TSV_BBVIA_30_FM1_1_BM1_1Xmin4 = 4.5*mm + TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Xmax4 = 5.25*mm - TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Ymin4 = -4.05*mm + 5*TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Ymax4 = 4.05*mm - 5*TSV_BBVIA_30_FM1_1_BM1_1R;
G4double TSV_BBVIA_30_FM1_1_BM1_1Xstep4 = (TSV_BBVIA_30_FM1_1_BM1_1Xmax4 - TSV_BBVIA_30_FM1_1_BM1_1Xmin4)/3;
G4double TSV_BBVIA_30_FM1_1_BM1_1Ystep4 = (TSV_BBVIA_30_FM1_1_BM1_1Ymax4 - TSV_BBVIA_30_FM1_1_BM1_1Ymin4)/32;
for (G4int i=0; i<132; i++){
    G4double x = TSV_BBVIA_30_FM1_1_BM1_1Xmin4 + TSV_BBVIA_30_FM1_1_BM1_1Xstep4*(i%4);
    G4double y = TSV_BBVIA_30_FM1_1_BM1_1Ymin4 + TSV_BBVIA_30_FM1_1_BM1_1Ystep4*(i/4);
    new G4PVPlacement(0,
                      G4ThreeVector(x, y, 0.5 * (FM1_1High + BM1_1High)),
                      logicTSV_BBVIA_30_FM1_1_BM1_1,
                      "TSV_BBVIA_30_FM1_1_BM1_1",
                      logicEnv,
                      false,
                      nTSV_BBVIA_30_FM1_1_BM1_1++,
                      true);
}


    

//FP1_1,type:Dielectric,Material:Polyimide
G4double FP1_1High = 2.4454*mm;
G4Box* solidFP1_1 =
    new G4Box("FP1_1",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(FP1_1High-FM1_1High));
G4LogicalVolume* logicFP1_1 =
    new G4LogicalVolume(solidFP1_1,
                        PolyimideResin,
                        "FP1_1");
G4ThreeVector posFP1_1 = G4ThreeVector(0,0,0.5*(FP1_1High+FM1_1High));
    new G4PVPlacement(0,
                      posFP1_1,
                      logicFP1_1,
                      "FP1_1",
                      logicEnv,
                      false,
                      0,
                      true);

//UBM_1,type:Metal,Material:Filled
G4double UBM_1High = 2.4494*mm;    
G4Box* solidUBM_1 =
    new G4Box("UBM_1",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UBM_1High-FP1_1High));
G4LogicalVolume* logicUBM_1 =
    new G4LogicalVolume(solidUBM_1,
                        Filled,
                        "UBM_1");
G4ThreeVector posUBM_1 = G4ThreeVector(0,0,0.5*(UBM_1High+FP1_1High));
    new G4PVPlacement(0,
                      posUBM_1,
                      logicUBM_1,
                      "UBM_1",
                      logicEnv,
                      false,
                      0,
                      true);

//焊盘
G4double TSVIPDPADR = 0.03*mm;
G4Tubs* solidTSVIPDPAD =
    new G4Tubs("TSVIPDPAD",
    0,
    TSVIPDPADR,
    0.5*(UBM_1High-FM1_1High),
    0.*deg,
    360.*deg);
G4LogicalVolume* logicTSVIPDPAD =
    new G4LogicalVolume(solidTSVIPDPAD,
                        copper,
                        "TSVIPDPAD");
G4int nTSVIPDPAD = 0;                        
for(G4double x = -4*mm; x<=4*mm; x+=0.5*mm){
    for(G4double y = 0.625*mm; y<=3.625*mm; y+=0.5*mm){
        G4double z = 0.5 * (FM1_1High + UBM_1High);
        if (x == 0.*mm&&(y==1.625*mm||y==1.125*mm)) continue;
        G4ThreeVector posTSVIPDPAD = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVIPDPAD,
                      logicTSVIPDPAD,
                      "TSVIPDPAD",
                      logicEnv,
                      false,
                      nTSVIPDPAD++,
                      true);
    }
}
for(G4double x = -4*mm; x<=4*mm; x+=0.5*mm){
    for(G4double y = -3.625*mm; y<=-0.125*mm; y+=0.5*mm){
        G4double z = 0.5 * (FM1_1High + UBM_1High);
        G4ThreeVector posTSVIPDPAD = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVIPDPAD,
                      logicTSVIPDPAD,
                      "TSVIPDPAD",
                      logicEnv,
                      false,
                      nTSVIPDPAD++,
                      true);
    }
}
for(G4double x = -5.125*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = 4.125*mm; y<=4.625*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_1High + UBM_1High);
        G4ThreeVector posTSVIPDPAD = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVIPDPAD,
                      logicTSVIPDPAD,
                      "TSVIPDPAD",
                      logicEnv,
                      false,
                      nTSVIPDPAD++,
                      true);
    }
}

for(G4double x = -5.125*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = -4.625*mm; y<=-4.125*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_1High + UBM_1High);
        G4ThreeVector posTSVIPDPAD = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVIPDPAD,
                      logicTSVIPDPAD,
                      "TSVIPDPAD",
                      logicEnv,
                      false,
                      nTSVIPDPAD++,
                      true);
    }
}

for(G4double x = -5.125*mm; x<=-4.625*mm; x+=0.25*mm){
    for(G4double y = -3.875*mm; y<=3.875*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_1High + UBM_1High);
        G4ThreeVector posTSVIPDPAD = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVIPDPAD,
                      logicTSVIPDPAD,
                      "TSVIPDPAD",
                      logicEnv,
                      false,
                      nTSVIPDPAD++,
                      true);
    }
}

for(G4double x = 4.625*mm; x<=5.125*mm; x+=0.25*mm){
    for(G4double y = -3.875*mm; y<=3.875*mm; y+=0.25*mm)    {
        G4double z = 0.5 * (FM1_1High + UBM_1High);
        G4ThreeVector posTSVIPDPAD = G4ThreeVector(x, y, z);
        new G4PVPlacement(0,
                      posTSVIPDPAD,
                      logicTSVIPDPAD,
                      "TSVIPDPAD",
                      logicEnv,
                      false,
                      nTSVIPDPAD++,
                      true);
    }
}

//UNNAMED_002,type:Metal,Material:FR_4
G4double UNNAMED_002High = 2.4574*mm;
G4Box* solidUNNAMED_002 =
    new G4Box("UNNAMED_002",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(UNNAMED_002High-UBM_1High));
G4LogicalVolume* logicUNNAMED_002 =
    new G4LogicalVolume(solidUNNAMED_002,
                        EpoxyResin,
                        "UNNAMED_002");
G4ThreeVector posUNNAMED_002 = G4ThreeVector(0,0,0.5*(UNNAMED_002High+UBM_1High));
    new G4PVPlacement(0,
                      posUNNAMED_002,
                      logicUNNAMED_002,
                      "UNNAMED_002",
                      logicEnv,
                      false,
                      0,
                      true);

//DIE_PAD,type:Metal,Material:Filled
G4double DIE_PADHigh = 2.4614*mm;    
G4Box* solidDIE_PAD =
    new G4Box("DIE_PAD",
    0.5*chip_sizeX,0.5*chip_sizeY,0.5*(DIE_PADHigh-UNNAMED_002High));
G4LogicalVolume* logicDIE_PAD =
    new G4LogicalVolume(solidDIE_PAD,
                        Filled,
                        "DIE_PAD");
G4ThreeVector posDIE_PAD = G4ThreeVector(0,0,0.5*(DIE_PADHigh+UNNAMED_002High));
    new G4PVPlacement(0,
                      posDIE_PAD,
                      logicDIE_PAD,
                      "DIE_PAD",
                      logicEnv,
                      false,
                      0,
                      true);           
fScoringVolume= logicEnv;

return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
