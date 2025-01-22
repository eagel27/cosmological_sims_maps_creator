QUERY_SUBHALO_DATA = """
    SELECT
        PROG.GalaxyID as ID,
        PROG.GroupNumber as GN,
        PROG.SubGroupNumber as SGN,
        PROG.Redshift,
        PROG.LastProgID as LastProgID,
        PROG.TopLeafID as TopLeafID,
        PROG.DescendantID as DescendantID,
        PROG.CentreOfPotential_x as xCM,
        PROG.CentreOfPotential_y as yCM,
        PROG.CentreOfPotential_z as zCM,
        PROG.Velocity_x as Vx,
        PROG.Velocity_y as Vy,
        PROG.Velocity_z as Vz,
        PROG.MassType_Star as MS,
        PROG.HalfMassRad_Star as HMRS
    FROM
        %s_SubHalo as PROG
    WHERE
        PROG.SnapNum=%s
        and PROG.MassType_Star > 1e10
    ORDER BY
        PROG.MassType_Star desc
"""


MPB_QUERY = ("""
    SELECT         
            REF.GalaxyID,       
            REF.MassType_Star,                  
            MPB.GalaxyID as MPB_GalaxyID,       
            MPB.SnapNum as MPB_SnapNum      
    FROM       
            %s_SubHalo as REF,       
            %s_SubHalo as MPB      
    WHERE       
            REF.SnapNum=%s       
            and REF.MassType_Star > 1e10         
            and MPB.GalaxyID between REF.GalaxyID and REF.TopLeafID        
    ORDER BY       
            REF.MassType_Star desc,       
            MPB.SnapNum desc    
""",
             2,
             ('GalaxyID', 'Mass',
              'MPB_GalaxyID', 'MPB_SnapNum'))


GALAXY_IDS_QUERY = ("""
            SELECT                           
                    REF.GalaxyID,               
                    REF.MassType_Star                      
            FROM               
                    %s_SubHalo as REF              
            WHERE               
                    REF.SnapNum=28      
                    and REF.MassType_Star > 1e9                   
            ORDER BY               
                    REF.MassType_Star desc            
""",
                    1,
                    ('GalaxyID', 'Mass'))


GALAXY_PARAMETERS = ("""
        SELECT                
             REF.GalaxyID,   
             REF.HalfMassRad_Star, 
             REF.MassType_Star,  
             REF.MassType_Gas,
             REF.MassType_DM,
             REF.MassType_BH,
             REF.GasSpin_x,
             REF.GasSpin_y,
             REF.GasSpin_z,
             REF.Stars_Spin_x, 
             REF.Stars_Spin_y, 
             REF.Stars_Spin_z, 
             REF.GroupNumber,
             REF.SubGroupNumber,
             SIZE.R_halfmass100,  
             SIZE.R_halfmass30                     
        FROM                
             %s_SubHalo as REF,       
             %s_Sizes as SIZE               
        WHERE                
             REF.SnapNum=%s                
             and REF.MassType_Star > 1e9          
             and REF.GalaxyID=SIZE.GalaxyID      
        ORDER BY                
             REF.MassType_Star desc
""",
                     2,
                     ('GalaxyID', 'Mass_Star', 'Mass_Gas', 'Mass_DM', 'Mass_BH', 'HalfMassRad_Star',
                      'R_halfmass100', 'R_halfmass30'),
                     "galaxy_params")

