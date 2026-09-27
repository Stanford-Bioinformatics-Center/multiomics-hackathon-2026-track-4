//Tissue
CREATE
  (human_muscle:Tissue
    {tissue_id: 'human_muscle', species: 'Human', name: 'Skeletal muscle'}),
  (human_adipose:Tissue
    {tissue_id: 'human_adipose', species: 'Human', name: 'Adipose'}),
  (human_plasma:Tissue
    {tissue_id: 'human_plasma', species: 'Human', name: 'Plasma'}),
  (human_whole_blood:Tissue
    {tissue_id: 'human_whole-blood', species: 'Human', name: 'Whole blood'}),
  (rat_ADRNL:Tissue
    {tissue_id: 'rat_ADRNL', species: 'Rat', name: 'Adrenal gland'}),
  (rat_BAT:Tissue
    {tissue_id: 'rat_BAT', species: 'Rat', name: 'Brown adipose'}),
  (rat_BLOOD:Tissue {tissue_id: 'rat_BLOOD', species: 'Rat', name: 'Blood'}),
  (rat_COLON:Tissue {tissue_id: 'rat_COLON', species: 'Rat', name: 'Colon'}),
  (rat_CORTEX:Tissue
    {tissue_id: 'rat_CORTEX', species: 'Rat', name: 'Cerebral cortex'}),
  (rat_HEART:Tissue {tissue_id: 'rat_HEART', species: 'Rat', name: 'Heart'}),
  (rat_HIPPOC:Tissue
    {tissue_id: 'rat_HIPPOC', species: 'Rat', name: 'Hippocampus'}),
  (rat_HYPOTH:Tissue
    {tissue_id: 'rat_HYPOTH', species: 'Rat', name: 'Hypothalamus'}),
  (rat_KIDNEY:Tissue {tissue_id: 'rat_KIDNEY', species: 'Rat', name: 'Kidney'}),
  (rat_LIVER:Tissue {tissue_id: 'rat_LIVER', species: 'Rat', name: 'Liver'}),
  (rat_LUNG:Tissue {tissue_id: 'rat_LUNG', species: 'Rat', name: 'Lung'}),
  (rat_OVARY:Tissue {tissue_id: 'rat_OVARY', species: 'Rat', name: 'Ovary'}),
  (rat_PLASMA:Tissue {tissue_id: 'rat_PLASMA', species: 'Rat', name: 'Plasma'}),
  (rat_SKM_GN:Tissue
    {tissue_id: 'rat_SKM-GN', species: 'Rat', name: 'Gastrocnemius muscle'}),
  (rat_SKM_VL:Tissue
    {tissue_id: 'rat_SKM-VL', species: 'Rat', name: 'Vastus lateralis muscle'}),
  (rat_SMLINT:Tissue
    {tissue_id: 'rat_SMLINT', species: 'Rat', name: 'Small intestine'}),
  (rat_SPLEEN:Tissue {tissue_id: 'rat_SPLEEN', species: 'Rat', name: 'Spleen'}),
  (rat_TESTES:Tissue {tissue_id: 'rat_TESTES', species: 'Rat', name: 'Testes'}),
  (rat_VENACV:Tissue
    {tissue_id: 'rat_VENACV', species: 'Rat', name: 'Vena cava'}),
  (rat_WAT_SC:Tissue
    {
      tissue_id: 'rat_WAT-SC',
      species: 'Rat',
      name: 'Subcutaneous white adipose'
    });

//Gene
CREATE
  CREATE (SLC7A5:Gene {symbol: 'SLC7A5','ncbi':8140,'ensembl':ENSG00000103257}),
  (SLC38A4:Gene {symbol: 'SLC38A4','ncbi':55089,'ensembl':ENSG00000139209}),
  (ALDH7A1:Gene {symbol: 'ALDH7A1','ncbi':501,'ensembl':ENSG00000164904}),
  (OARD1:Gene {symbol: 'OARD1','ncbi':221443,'ensembl':ENSG00000124596}),
  (PNP:Gene {symbol: 'PNP','ncbi':4860,'ensembl':ENSG00000198805}),
  (CHKB:Gene {symbol: 'CHKB','ncbi':1120,'ensembl':ENSG00000100288}),
  (BCAT2:Gene {symbol: 'BCAT2','ncbi':587,'ensembl':ENSG00000105552}),
  (AADAT:Gene {symbol: 'AADAT','ncbi':51166,'ensembl':ENSG00000109576}),
  (GPCPD1:Gene {symbol: 'GPCPD1','ncbi':56261,'ensembl':ENSG00000125772}),
  (ETNK1:Gene {symbol: 'ETNK1','ncbi':55500,'ensembl':ENSG00000139163}),
  (Slc7a5:Gene {symbol: 'Slc7a5','ensembl':ENSRNOG00000018824}),
  (Slc38a4:Gene {symbol: 'Slc38a4','ensembl':ENSRNOG00000006653}),
  (Aldh7a1:Gene {symbol: 'Aldh7a1','ensembl':ENSRNOG00000014645}),
  (Oard1:Gene {symbol: 'Oard1','ensembl':ENSRNOG00000012679}),
  (Pnp:Gene {symbol: 'Pnp','ensembl':ENSRNOG00000009982}),
  (Chkb:Gene {symbol: 'Chkb','ensembl':ENSRNOG00000011404}),
  (Bcat2:Gene {symbol: 'Bcat2','ensembl':ENSRNOG00000020956}),
  (Aadat:Gene {symbol: 'Aadat','ensembl':ENSRNOG00000011861}),
  (Gpcpd1:Gene {symbol: 'Gpcpd1','ensembl':ENSRNOG00000053201}),
  (Etnk1:Gene {symbol: 'Etnk1','ensembl':ENSRNOG00000014856});

//Protein

//Metabolite
CREATE
  (Xanthosine:Metabolite
    {
      metabolite_id: 'REFMET:Xanthosine',
      refmet_name: 'Xanthosine',
      metabolite_name: 'Xanthosine'
    });
CREATE
  (Hypoxanthine:Metabolite
    {
      metabolite_id: 'REFMET:Hypoxanthine',
      refmet_name: 'Hypoxanthine',
      metabolite_name: 'Hypoxanthine'
    });
CREATE
  (Malic_acid:Metabolite
    {
      metabolite_id: 'REFMET:Malic acid',
      refmet_name: 'Malic acid',
      metabolite_name: 'Malic acid'
    });
CREATE
  (Lactic_acid:Metabolite
    {
      metabolite_id: 'REFMET:Lactic acid',
      refmet_name: 'Lactic acid',
      metabolite_name: 'Lactic acid'
    });
CREATE
  (Succinic_acid:Metabolite
    {
      metabolite_id: 'REFMET:Succinic acid',
      refmet_name: 'Succinic acid',
      metabolite_name: 'Succinic acid'
    });
CREATE
  (Xanthine:Metabolite
    {
      metabolite_id: 'REFMET:Xanthine',
      refmet_name: 'Xanthine',
      metabolite_name: 'Xanthine'
    });
CREATE
  (Alanine:Metabolite
    {
      metabolite_id: 'REFMET:Alanine',
      refmet_name: 'Alanine',
      metabolite_name: 'Alanine'
    });
CREATE
  (Inosine:Metabolite
    {
      metabolite_id: 'REFMET:Inosine',
      refmet_name: 'Inosine',
      metabolite_name: 'Inosine'
    });
CREATE
  (Pantothenic_acid:Metabolite
    {
      metabolite_id: 'REFMET:Pantothenic acid',
      refmet_name: 'Pantothenic acid',
      metabolite_name: 'Pantothenic acid'
    });
CREATE
  (Pyruvic_acid:Metabolite
    {
      metabolite_id: 'REFMET:Pyruvic acid',
      refmet_name: 'Pyruvic acid',
      metabolite_name: 'Pyruvic acid'
    });