import numpy  as np
import pandas as pd

def _guess_pairs14(bonds):
    pairs14 = []
    
    #if len(self.dsorted) == 0:
    #    return pairs14

    for b in bonds:
        for i in [0, 1]:
            idx1 = b[i]
            bA1l = bonds[:,0] == idx1
            bA1r = bonds[:,1] == idx1

            ids2 = np.concatenate([bonds[:,1][bA1l], bonds[:,0][bA1r]])
            if len(ids2) == 0: continue
            
            bA2l = np.isin(bonds[:,0], ids2)
            bA2r = np.isin(bonds[:,1], ids2)

            ids3 = np.concatenate([bonds[:,1][bA2l], bonds[:,0][bA2r]])
            ids3 = np.delete(ids3, np.isin(ids3, idx1))
            if len(ids3) == 0: continue
            
            bA3l = np.isin(bonds[:,0], ids3)
            bA3r = np.isin(bonds[:,1], ids3)

            ids4 = np.concatenate([bonds[:,1][bA3l], bonds[:,0][bA3r]])
            ids4 = np.delete(ids4, np.isin(ids4, ids2))
            ids4 = np.delete(ids4, np.isin(ids4, ids3)) #cholesterol pentagon
            if len(ids4) == 0: continue

            for idx4 in ids4:
                if idx1 > idx4:
                    result = [idx4, idx1]

                else:
                    result = [idx1, idx4]

                if result not in pairs14:
                    pairs14.append(result)

    df = pd.DataFrame(pairs14)
    return df.sort_values(by=[0,1]).to_numpy()

