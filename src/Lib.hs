module Lib where

import Data.Maybe (fromMaybe)
import Numeric.GSL.Integration
import Test.LeanCheck
import Text.PrettyPrint

-- Doc Utility
----------------------------------------------------

eqDoc :: String -> String -> Doc
eqDoc name quantity = hsep $ map text
    [ name, "=", quantity ]
    
capIndent :: String -> [Doc] -> Doc
capIndent cap ds = hang (text cap) 4 $ vcat ds
    
-- Cylinder Params
----------------------------------------------------
data CylinderParams = CylinderParams
    { parL, parR, parH0, parHL :: Double }
    deriving Show
    
-- Preconditions on our parameters
--      0.0 < L < Infinity
--      0.0 < r < Infinity
--      h_0, h_L \in [0, 2r]
------------------------------------------
preCond_L :: CylinderParams -> Bool
preCond_L p =
        0.0 < len
    &&  (not $ isInfinite len)
    where
        len = parL p

preCond_r :: CylinderParams -> Bool
preCond_r p =
        0.0 < r
    &&  (not $ isInfinite r)
    where
        r = parR p

preCond_hs :: CylinderParams -> Bool
preCond_hs p =
        (0.0 <= h_0 && h_0 <= 2*r)
    &&  (0.0 <= h_L && h_L <= 2*r)
    &&  h_0 /= h_L
    where
        r = parR p
        h_0 = parH0 p
        h_L = parHL p

preCond_Params :: CylinderParams -> Bool
preCond_Params p =
        preCond_L p 
    &&  preCond_r p 
    &&  preCond_hs p
    
-- Utility functions
------------------------------------------
withHeights :: Double -> Double -> CylinderParams
withHeights h_0 h_L = CylinderParams
    { parL = 1, parR = 1, parH0 = h_0, parHL = h_L }
    
flipHeights :: CylinderParams -> CylinderParams
flipHeights p = CylinderParams 
    {   parL = parL p, parR = parR p, 
        parH0 = parHL p, parHL = parH0 p }
        
paramsDoc :: CylinderParams -> Doc
paramsDoc p = capIndent
    "Cylinder parameters:"
    [ eqDoc "L  " $ show $ parL p
    , eqDoc "r  " $ show $ parR p
    , eqDoc "h_0" $ show $ parH0 p
    , eqDoc "h_L" $ show $ parHL p
    ]
        
printParams :: CylinderParams -> IO ()
printParams p = do
    putStrLn $ render $ paramsDoc p

-- Functions
----------------------------------------------------
    
sigma :: Double -> Double -> Double
sigma r x = sqrt $ 2*r*x - x^2

heightFunction :: Double -> Double -> Double -> 
    Double -> Double
heightFunction h_0 h_L len x =
    ((len - x) * h_0 + x * h_L) / len

areaFunction :: Double -> Double -> Double
areaFunction r h = r^2 * 
    acos (1 - h / r) - (r - h) * sigma r h

volumeFormula :: CylinderParams -> Double
volumeFormula p = (r^2/c_h) * ((rho h_0) - (rho h_L) + 
    ((cr h_L)*(s h_L) - (cr h_0)*(s h_0))/3)
    where
        len = parL p
        r = parR p
        h_0 = parH0 p
        h_L = parHL p
        
        c_h = (h_L - h_0) / len
        
        rho :: Double -> Double
        rho x = (r - x) * acos((r - x)/r)
        
        cr :: Double -> Double
        cr x = 2 + ((r - x)/r)^2
        
        s :: Double -> Double
        s = sigma r

-- Comparison Formula <-> Numeric
----------------------------------------------------
type GSLResult = (Double, Double)   -- (Result, Error)

gslVolume :: CylinderParams -> GSLResult
gslVolume p = integrateQNG prec (a . h) 0 len
    where
        prec = 1e-9
        len = parL p
        r = parR p
        h_0 = parH0 p
        h_L = parHL p
        
        h :: Double -> Double
        h = heightFunction h_0 h_L len
        
        a :: Double -> Double
        a = areaFunction r

gslResultsDoc :: CylinderParams -> Double -> GSLResult -> Doc
gslResultsDoc p formula (estimate, error) = vcat
    [ pDoc, fDoc, eDoc, dDoc ]
    where
        pDoc = paramsDoc p

        fDoc = capIndent
            "Formula result:"
            [ eqDoc "V" $ show formula ]

        eDoc = capIndent
            "Numeric Integration result:"
            [ eqDoc "V" $ show estimate
            , eqDoc "Error" $ show error
            ]

        dDoc = capIndent
            "Difference Formula - Best Estimate:"
            [ eqDoc "V_f - V_e" $ show $ formula - estimate ]

compareFormulaToGSL :: CylinderParams -> IO ()
compareFormulaToGSL p = 
    putStrLn $ render $ gslResultsDoc p formula gsl
    where
        formula = volumeFormula p
        gsl = gslVolume p

-- LeanCheck
----------------------------------------------------
instance Listable CylinderParams where
    tiers = (cons4 CylinderParams) `suchThat` preCond_Params

withinTolerance :: Double -> Double -> Double -> Bool
withinTolerance delta x y = abs (x - y) < delta

prop_VolumeWithinToleranceGSL :: Double -> 
    CylinderParams -> Bool
prop_VolumeWithinToleranceGSL delta p =
    withinTolerance delta formula estimate
    where
        formula = volumeFormula p
        (estimate, _) = gslVolume p

test0 = checkFor 10000 $ prop_VolumeWithinToleranceGSL 1e-8
test1 = checkFor 100000 $ prop_VolumeWithinToleranceGSL 1e-8
