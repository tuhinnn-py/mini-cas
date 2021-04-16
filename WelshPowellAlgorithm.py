def use_welsh_powell(adj):
    __colors, N = 0, len(adj)
    __valence = [[0, i] for i in range(N)]
    __color_dict = [None for i in range(N)]

    for i in range(N):
        for j in range(N):
            if not adj[i][j] == None:
                __valence[i][0] += 1
                __valence[j][0] += 1

    __valence = sorted(__valence, reverse = True, key = lambda x : x[0])
    
    while __valence:
        __node = __valence.pop(0)[1]
        if not __color_dict[__node] == None:
            continue
        
        __colors += 1
        __color_dict[__node] = __colors
        __colored = set([__node])

        for i in __valence:
            i, __safe = i[1], True
            
            for j in range(N):
                if not(adj[i][j] == None) and j in __colored:
                    __safe = False
                    break
                
            if __safe:
                __color_dict[i] = __colors
                __colored.add(i)
        
    return __colors, __color_dict
