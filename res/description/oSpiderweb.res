CONTAINER oSpiderweb
{
    NAME oSpiderweb;
    INCLUDE Obase;

    GROUP ID_OBJECTPROPERTIES
    {
        BOOL USE_RAYCAST {}
        REAL RADIUS_WEB {MIN 0.0;}
        LONG ANCHOR_POINTS_CNT {MIN 3; MAX 9;}
        IN_EXCLUDE OBJECT_LIST {}

        LONG SUBDIVISION {MIN 1;}
        LONG SEED_ID {MIN 1;}

        GROUP ID_GRP_RADIAL 
        {
            DEFAULT 1; 
            LONG RADII_CNT {MIN 1;}
            REAL RND_CENTER_DISTANCE {UNIT PERCENT;MIN 0.0; MAX 100.0;}
        }
        GROUP ID_GRP_SPIRAL
        {
            DEFAULT 1 ;
            GROUP ID_GRP_OUT_SPIRAL
            {
                DEFAULT 1; 
                REAL START_DISTANCE{MIN 0.0;}
                LONG OS_REVOLUTION_CNT {MIN 0;}
                LONG OS_OFFSET_START {MIN 0;}
                REAL OS_SPACE {MIN 0.0;}
                GROUP 
                {
                    COLUMNS 2;
                    REAL OS_TURBULENCE {MIN 0.0;}
                    REAL OS_LENGTH_WEIGHT {UNIT PERCENT;MIN 0.0; MAX 100.0;}
                }
            }

            GROUP ID_GRP_CENTER_SPIRAL
            {
                DEFAULT 1; 
                BOOL CREATE_CENTRAL_SPIRAL{MIN 0;}
                LONG CS_REVOLUTION_CNT {MIN 0;}
                LONG CS_OFFSET_START {MIN 0;}
                REAL CS_SPACE {MIN 0;}
                GROUP
                {
                    COLUMNS 2;
                    REAL CS_TURBULENCE {MIN 0.0;}
                    REAL CS_LENGTH_WEIGHT {UNIT PERCENT;MIN 0.0; MAX 100.0;}
                }
            }


        }
        GROUP ID_GRP_DISPLAY
        {
            DEFAULT 1; 
            BOOL DRAWPHERE{}
            BOOL DRAWPLANE{}

        }


    }
}
