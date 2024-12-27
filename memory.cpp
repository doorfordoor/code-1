//引线
G4double BM1_1linewidth = 0.02*mm;
G4double BM1_1linelongth = 9.7*mm;  
G4double BM1_1linehigh = 0.004*mm;
G4Box* solidBM1_1line =
    new G4Box("BM1_1line",
    0.5*BM1_1linewidth, 0.5*BM1_1linelongth,0.5*BM1_1linehigh);  
G4LogicalVolume* logicBM1_1line =
    new G4LogicalVolume(solidBM1_1line,
                        copper,
                        "BM1_1line");
G4int nBM1_1line = 0;
for(G4double x = -5.285*mm+0.5*BM1_1linewidth; x<=1.142*mm-0.5*BM1_1linewidth; x+=0.05*mm){
    G4ThreeVector posBM1_1line = G4ThreeVector(x,0,0);
    new G4PVPlacement(0,
                      posBM1_1line,
                      logicBM1_1line,
                      "BM1_1line",
                      logicBM1_1,
                      false,
                      nBM1_1line++,
                      true);
}                    
G4double BM1_1linelongth1 = 4.835*mm;
G4Box* solidBM1_1line1 =
    new G4Box("BM1_1line1",
    0.5*BM1_1linewidth, 0.5*BM1_1linelongth1,0.5*BM1_1linehigh);
G4LogicalVolume* logicBM1_1line1 =
    new G4LogicalVolume(solidBM1_1line1,
                        copper,
                        "BM1_1line1");
G4double nBM1_1line1 = 0;
for(G4double x = 1.142*mm+0.5*BM1_1linewidth; x<=3.5216*mm-0.5*BM1_1linewidth; x+=0.05*mm){
    G4ThreeVector posBM1_1line1 = G4ThreeVector(x,-2.4025*mm,0);  
    new G4PVPlacement(0,
                      posBM1_1line1,
                      logicBM1_1line1,
                      "BM1_1line1",
                      logicBM1_1,
                      false,
                      nBM1_1line1++,
                      true);                
}
G4double BM1_1linelongth2 = 2.025*mm;    
G4Box* solidBM1_1line2 =
    new G4Box("BM1_1line2",
    0.5*BM1_1linewidth, 0.5*BM1_1linelongth2,0.5*BM1_1linehigh);
G4LogicalVolume* logicBM1_1line2 =
    new G4LogicalVolume(solidBM1_1line2,
                        copper,
                        "BM1_1line2");
G4double nBM1_1line2 = 0;
for(G4double x = 1.142*mm+0.5*BM1_1linewidth; x<=4.14*mm-0.5*BM1_1linewidth; x+=0.05*mm){
    G4ThreeVector posBM1_1line2 = G4ThreeVector(x,3.8375*mm,0); 
    new G4PVPlacement(0,
                      posBM1_1line2,
                      logicBM1_1line2,
                      "BM1_1line2",
                      logicBM1_1,
                      false,
                      nBM1_1line2++,
                      true);
}
G4double BM1_1linelongth3 = 3.79*mm; 
G4Box* solidBM1_1line3 =
    new G4Box("BM1_1line3",
    0.5*BM1_1linewidth, 0.5*BM1_1linelongth3,0.5*BM1_1linehigh);
G4LogicalVolume* logicBM1_1line3 =
    new G4LogicalVolume(solidBM1_1line3,
                        copper,
                        "BM1_1line3");
G4double nBM1_1line3 = 0;
for(G4double x = 4.155*mm+0.5*BM1_1linewidth; x<=5.285*mm-0.5*BM1_1linewidth; x+=0.05*mm){
    G4ThreeVector posBM1_1line3 = G4ThreeVector(x,2.955*mm,0);
    new G4PVPlacement(0,
                      posBM1_1line3,
                      logicBM1_1line3,
                      "BM1_1line3",
                      logicBM1_1,
                      false,
                      nBM1_1line3++,
                      true);
}
G4Box *solidBM1_1line4 =
    new G4Box("BM1_1line4",
    0.5*0.875*mm, 0.5*0.54083*mm,0.5*BM1_1linehigh);
G4LogicalVolume* logicBM1_1line4 =
    new G4LogicalVolume(solidBM1_1line4,
                        copper,
                        "BM1_1line4");
new G4PVPlacement(0,
                      G4ThreeVector(4.8425*mm,0.7617*mm,0),
                      logicBM1_1line4,
                      "BM1_1line4",
                      logicBM1_1,
                      false,
                      0,
                      true);   
G4Box *solidBM1_1line5 =
    new G4Box("BM1_1line5",
    0.5*0.875*mm, 0.5*0.4405*mm,0.5*BM1_1linehigh);
G4LogicalVolume* logicBM1_1line5 =
    new G4LogicalVolume(solidBM1_1line5,
                        copper,
                        "BM1_1line5");
new G4PVPlacement(0,
                      G4ThreeVector(4.8425*mm,0.25105*mm,0),
                      logicBM1_1line5,
                      "BM1_1line5",
                      logicBM1_1,
                      false,
                      0,
                      true);
G4Box *solidBM1_1line6 =
    new G4Box("BM1_1line6",
    0.5*1.42*mm, 0.5*0.8755*mm,0.5*BM1_1linehigh);     
G4LogicalVolume* logicBM1_1line6 =
    new G4LogicalVolume(solidBM1_1line6,
                        copper,
                        "BM1_1line6");
new G4PVPlacement(0,
                      G4ThreeVector(4.23*mm,-4.41225*mm,0),
                      logicBM1_1line6,
                      "BM1_1line6",
                      logicBM1_1,
                      false,
                      0,
                      true);                 



//FM1_1line,copper,导线
G4double FM1_1linelongth = 9.7*mm;
G4double FM1_1linewidth = 0.02*mm;
G4double FM1_1linehigh = 0.004*mm;
G4Box *solidFM1_1line =
    new G4Box("FM1_1line",
    0.5*FM1_1linewidth,0.5*FM1_1linelongth,0.5*FM1_1linehigh);
G4LogicalVolume* logicFM1_1line =
    new G4LogicalVolume(solidFM1_1line,
                        copper,
                        "FM1_1line");
G4int nFM1_1line = 0;
for(G4double x = -5.285*mm + 0.5*FM1_1linewidth; x < 5.285*mm - 0.5*FM1_1linewidth; x += 0.05*mm){
    G4ThreeVector posFM1_1line = G4ThreeVector(x,0,0);
    new G4PVPlacement(0,
                      posFM1_1line,
                      logicFM1_1line,
                      "FM1_1line",
                      logicFM1_1,
                      false,
                      nFM1_1line++,
                      true);
}

//FM1_1DDR3,引脚
G4double FM1_1DDR3width = 0.047*mm;
G4double FM1_1DDR3high = 0.059*mm;
G4Box *solidFM1_1DDR3 =
    new G4Box("FM1_1DDR3",
    0.5*FM1_1DDR3width,0.5*FM1_1DDR3width,0.5*FM1_1linehigh);
G4LogicalVolume* logicFM1_1DDR3 =
    new G4LogicalVolume(solidFM1_1DDR3,
                        copper,
                        "FM1_1DDR3Pin");
G4double FM1_1DDR3xmax1 = 4.26013*mm;  
G4double FM1_1DDR3xmin1 = 2.33726*mm;   
G4double FM1_1y = 0.1543*mm;    
G4double FM1_1DDR3xstep1 = (FM1_1DDR3xmax1 - FM1_1DDR3xmin1)/29;
for(G4int i = 0; i < 30; i++){
    G4double x = FM1_1DDR3xmin1 + FM1_1DDR3xstep1*i;
    G4ThreeVector posFM1_1DDR3 = G4ThreeVector(x,FM1_1y,0);
    new G4PVPlacement(0,
                      posFM1_1DDR3,
                      logicFM1_1DDR3,
                      "FM1_1DDR3",
                      logicFM1_1,
                      false,
                      i,
                      true);      
}      
G4double FM1_1DDR3xmax2 = 2.1068*mm;
G4double FM1_1DDR3xmin2 = -0.43972*mm;
G4double FM1_1DDR3xstep2 = (FM1_1DDR3xmax2 - FM1_1DDR3xmin2)/36;
for(G4int i = 0; i < 37; i++){
    G4double x = FM1_1DDR3xmin2 + FM1_1DDR3xstep2*i;
    G4ThreeVector posFM1_1DDR3 = G4ThreeVector(x,FM1_1y,0);
    new G4PVPlacement(0,
                      posFM1_1DDR3,
                      logicFM1_1DDR3,
                      "FM1_1DDR3",
                      logicFM1_1,
                      false,
                      i+30,
                      true);
}
G4double FM1_1DDR3xmax3 = -0.6371*mm;
G4double FM1_1DDR3xmin3 = -2.01868*mm;
G4double FM1_1DDR3xstep3 = (FM1_1DDR3xmax3 - FM1_1DDR3xmin3)/21;
for(G4int i = 0; i < 22; i++){
    G4double x = FM1_1DDR3xmin3 + FM1_1DDR3xstep3*i;
    G4ThreeVector posFM1_1DDR3 = G4ThreeVector(x,FM1_1y,0);
    new G4PVPlacement(0,
                      posFM1_1DDR3,
                      logicFM1_1DDR3,
                      "FM1_1DDR3",
                      logicFM1_1,
                      false,
                      i+67,
                      true);
}
G4double FM1_1DDR3xmax4 = -2.35222*mm;
G4double FM1_1DDR3xmin4 = -4.26013*mm;
G4double FM1_1DDR3xstep4 = (FM1_1DDR3xmax4 - FM1_1DDR3xmin4)/28;
for(G4int i = 0; i < 29; i++){
    G4double x = FM1_1DDR3xmin4 + FM1_1DDR3xstep4*i;
    G4ThreeVector posFM1_1DDR3 = G4ThreeVector(x,FM1_1y,0);
    new G4PVPlacement(0,
                      posFM1_1DDR3,
                      logicFM1_1DDR3,
                      "FM1_1DDR3",
                      logicFM1_1,
                      false,
                      i+89,
                      true);
}

