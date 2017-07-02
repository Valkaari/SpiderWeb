import c4d
#Welcome to the world of Python

from c4d.utils import GeRayCollider
from c4d import utils
import random

##todo :
## add random to support thread
## random center for radial thread 
## 
## 
##
##  radiicnt = 26, amin=0, amax=2, bmin=3, bmax = 4, extSDist = None,nbTour = 50)


##
## @brief      Class Root for spider web generator
##
class SpiderGenerator:
    def __repr__(self):
        return str(self)

    def __init__(self, hub= None, subdiv = 8, seed  = 123456,radii = 26):
        
        if hub is None:
            return
        else:
            self.hub = hub
            self.hubMg = hub.GetMg()
            self.hubPos = self.hubMg.off
            self.subdiv = subdiv
            self.seed = seed
            self.radii = radii
            random.seed(seed)
            

    ##
    ## @brief      anchors can be a vector or a list of vectors. This function get the length of the vector or the first element in the list.
    ##
    ## @param      self  
    ## @param      a     the vector or the list to calculate the length from
    ##
    ## @return     
    ##
    def _getLengthVector(self,a):
        if type(a[1]) is list:
            return a[1][0].GetLength()
        else:
            return a[1].GetLength()
    ##
    ## @brief      anchors can be vector or list of vector. This return the vector or the first element of the list if it's a list
    ##
    ## @param      self  
    ## @param      v     the vector or the list of vector
    ##
    ## @return     
    ##
    def _getVector(self, v):
        if type(v) is list:
            return v[0]
        else:
            return v
        return None
    ##
    ## @brief      compute points positions using _getPointPositionAt() 
    ##
    ## @param      self  
    ## @param      p1    point start position vector
    ## @param      p2    point end position vector
    ## @param      hub   position of the hub
    ## @param      pcnt  number of point to compute
    ##
    ## @return     return array of points
    ##
    def _getPointsPosition(self,p1, p2 , hubPos, pcnt):
        #compute points on spline moving toward the hub
        points  = []
        for cpt in xrange(pcnt+1):
            t = float(cpt) / float(pcnt)
            points.append( self._getPointPositionAt(p1,p2,hubPos,t) )
        return points
    ##
    ## @brief      compute bezier curve 3 points
    ##
    ## @param      self  
    ## @param      p1    point start of the curve
    ## @param      p2    Point end of the curve
    ## @param      hub   position of the center of the spider web (were the bezier curve have to go)
    ## @param      t     position of the point in the curve 0.0 < t < 1.0
    ##
    ## @return     return the point position
    ##
    def _getPointPositionAt(self,p1,p2,hubPos,t):
        #calcul la position d'un points sur une spline bézier 
        p1 +=hubPos
        p2 +=hubPos
        mp =  (p2 + p1 ) * 0.5
        direction = hubPos - mp
        direction.Normalize()
        direction *= (p2-p1).GetLength() * 0.2
        mp += direction 
        
        #calcul des positions 
        m = pow(1-t,2) * p1 + 2*(1-t)*t*mp + pow(t,2)*p2
        return m - hubPos
    

##
## @brief      Main class for Threads of the spider web
##
class Thread (SpiderGenerator,object):
    ##
    ##
    ##
    def __init__(self, **kwargs):
        super (Thread,self).__init__(**kwargs)
##
## @brief      main class to compute and display the web
##
class ComputeAndDisplay(SpiderGenerator,object):
    def __init__(self, radius = 400, desiredAnchorsCnt = 3 ,  **kwargs):
        self.radius = radius
        self.desiredAnchorsCnt = desiredAnchorsCnt
        super (ComputeAndDisplay,self).__init__(**kwargs)

   
##
## @brief     SimpleThread is a thread composed by only two points, can be straight or curved.
##
class SimpleThread (Thread, object):
    ##
    ## @brief      constructor
    ##
    ## @param      self    
    ## @param      kwargs  arguments to call the parent
    ##
    def __init__(self, startPoint = None, endPoint = None, **kwargs):
        self.startPoint = startPoint
        self.endPoint = endPoint
        self.points = []
        super(SimpleThread, self).__init__(**kwargs)
    ##
    ## @brief      get points
    ##
    ## @param      self  
    ##
    ## @return     return list of  points
    ##
    def _getPoints(self):
        if len(self.points) > 0:
            return self.points
        else: 
            return None
    ##
    ## @brief      set point value
    ##
    ## @param      self       
    ## @param      startPoint  start point coordinate
    ## @param      endPoint    end point coordinate
    ##
    ## @return     
    ##
    def _setPoints(self,startPoint=None,endPoint=None):
        if startPoint is not None:
            self.startPoint = startPoint
        if endPoint is not None:
            self.endPoint  = endPoint
    ##
    ## @brief      compute intermediates points
    ##  
    ## @param      self  
    ##
    ## @return     
    ##
    def _computePoints(self, type =  c4d.SPLINETYPE_LINEAR):
        if type == c4d.SPLINETYPE_LINEAR:
            mp = (self.startPoint + self.endPoint) * 0.5
            self.points.append(self.startPoint)
            self.points.extend(self._getPointsPosition(self.startPoint,self.endPoint, mp, self.subdiv)[1:-1])
            self.points.append(self.endPoint)
            return True
        elif type == c4d.SPLINETYPE_BEZIER:

            self.points.append(self.startPoint)
            self.points.extend(self._getPointsPosition(self.startPoint,self.endPoint, self.hubPos,self.subdiv)[1:-1])
            self.points.append(self.endPoint)
            return True
        else:
            return False
        return True
##
## @brief      multiple thread is composed by multiple point or anchor points. (frame Thread, spiral)
##
class MultipleThread(Thread,object):
    def __init__(self,**kwargs):
        super(MultipleThread, self).__init__(**kwargs)


##
## @brief      Support Thread are in the corner of the spider web.
##
class SupportThread(SimpleThread,object):
    ##
    ## @brief      Init clas and call parent
    ##
    ## @param      self   
    ## @param      kwargs  used to call parent and it's constructor
    ##
    def __init__(self,**kwargs):
        self.done = False
        super(SupportThread, self).__init__(**kwargs)
    ##
    ## @brief      print the value of start and end point (from parent class)
    ##
    ## @param      self  
    ##
    ## @return     a string of start and end point coordinate and if this corner is done or not
    ##
    def __str__(self):
        return    '\r\n' + str(self.startPoint) + ',' + str(self.endPoint) + ',' + str(self.done)
    ##
    ## @brief      return the value of 'self.done'
    ##
    ## @param      self  
    ##
    ## @return     true or false 
    ##
    def _isDone(self):
        return self.done
    ##
    ## @brief      set the support thread as done
    ##
    ## @param      self  
    ##
    ## @return     
    ##
    def _setDone(self):
        self.done = True

    
##
## @brief      class for RadialThread More a Type
##
class RadialThread(SimpleThread,object):
    ##
    ## @brief      init radial thread
    ##
    ## @param      self     
    ## @param      ray_dir  ray direction
    ## @param      length   length of the ray
    ## @param      kwargs    call parent constructor
    ##
    def __init__(self, ray_dir = None, length = None,**kwargs):
        self.length = length
        self.ray_dir = ray_dir
        super(RadialThread, self).__init__(**kwargs)
    ##
    ## @brief      return the length of the thread
    ##
    ## @param      self  
    ##
    ## @return     the length of the thread
    ##
    def _getLength(self):
        return self.length
    ##
    ## @brief      set the length of the thread
    ##
    ## @param      self  
    ## @param      a     length of the thread
    ##
    ## @return     
    ##
    def _setLength(self,a):
        self.length = a
    ##
    ## @brief      return the point position at t with 0 < t <= 1
    ##
    ## @param      self  
    ## @param      t     point position 0 < t <= 1
    ##
    ## @return     { description_of_the_return_value }
    ##
    def _getPoint(self,t):

        return self.startPoint+ ((self.endPoint-self.startPoint) * t)
        
##
## @brief      The frame of the spider web
##
class FrameThread(MultipleThread,object):
    def __init__(self,**kwargs):
        self.anchors = []

        #call parent
        super(FrameThread,self).__init__(**kwargs)
    ##
    ## @brief      define anchors points and will calculate intermediate points at subdivision self.subdiv
    ##
    ## @param      self    
    ## @param      points  list of points
    ##               
    ## @return     True or False
    ##
    def _setPoints(self, points):
        if type(points) is list:
            self.anchors = list(points)
            self.points = []
            self._computePoints()
        else:
            return False
        return True

    ##
    ## @brief      Compute Points and intermediates points for the frame thread
    ##
    ## @param      self  
    ##
    ## @return     
    ##
    def _computePoints(self):
        for i, point in enumerate(self.anchors):
            apoint = self._getVector(point)
            npoint =  self._getVector(self.anchors[(i+1) % len(self.anchors)])
            self.points.append(apoint)
            self.points.extend(self._getPointsPosition(apoint,npoint,self.hubPos,self.subdiv))


    ##
    ## @brief      return the points list
    ##
    ## @param      self  
    ##
    ## @return     the points list
    ##
    def _getPoints(self):
        if len(self.points)>0:
            return self.points
        else:
            return None
        return None
    ##
    ## @brief      get the number of anchors
    ##
    ## @param      self  
    ##
    ## @return     the size of the array
    ##
    def _getAnchorCount(self):
        return len(self.anchors)
    ##
    ## @brief      return the number of composed anchors
    ##
    ## @param      self  
    ##
    ## @return     return the size of the array of the array created with only composed anchors
    ##
    def _getComposedAnchorCount(self):
        return len([t for t in self.anchors if type(t) is list])

    ##
    ## @brief      return only the anchors even if they are composed anchors.
    ##
    ## @param      self  
    ##
    ## @return     return a list of vectors (anchors)
    ##
    def _getAnchors(self):
        return [self._getVector(t) for t in self.anchors]
    ##
    ## @brief      return the composed anchors
    ##
    ## @param      self  
    ##
    ## @return     return a list of composed anchors with the sub points.
    ##
    def _getComposedAnchors(self):
        return [t for t in self.anchors if type(t) is list]

    def _getPointCount(self):
        return len(self.points)

##
## @brief      store the spiral part of the spiderweb
##
class SpiralThread(MultipleThread,object):
    ##
    ## @brief      The spiral is made of multiple segments of spine with xx points. 
    ##             init arrays to store point and number of point per segments.
    ##  
    ## @param      self    
    ## @param      kwargs  
    ##  
    def __init__(self,**kwargs):
        self.points = []
        self.segments = []

        #call parent
        super(SpiralThread,self).__init__(**kwargs)

    
    
##
## @brief      does the magic
##
class ComputeSpiderWeb(ComputeAndDisplay, object):
    ##
    ## @brief      Init values
    ##
    ## @param      self           
    ## @param      rotCnt         rotation number for the spiral thread
    ## @param      amin           a minimum value for spiral thread
    ## @param      amax           a maximum value for spiral thread
    ## @param      bmin           b minimum value for spiral thread
    ## @param      bmax           b maximum value for spiral thread
    ## @param      startDistance  start distance from center, if None it will be calculated
    ## @param      kwargs         other arguments for parent class
    ##
    def __init__(self, rotCnt = 50 , amin=0.0, amax=2.0, bmin=0.0, bmax = 4.0, startDistance = None, lengthWeight = 1.0, rndHubWeight = 1.0,offsetExtSpiral=0, 
                createCS = False, csRotCnt = 5 ,csamax = 2.0 , csbmax = 4.0, csOffSet =0 , cslgthWei = 1.0 ,**kwargs):
        
        super(ComputeSpiderWeb,self).__init__(**kwargs)
        #self.hairTag = hairTag
        self.objList =  []
        self.rcObjList = []
        ## new class
        self.frameThread   = FrameThread  (hub = self.hub, subdiv = self.subdiv)
        #self.supportThread = SupportThread(hub = self.hub, subdiv = self.subdiv, seed = self.seed)
        self.supportThreads =  []
        self.radialThreads  = []
        self.spiralThreads  = SpiralThread (hub = self.hub, subdiv = self.subdiv, seed= self.seed)
        self.rotCnt = rotCnt
        self.amin = amin
        self.amax = amax
        self.bmin = bmin
        self.bmax =bmax
        self.startDistance = startDistance
        self.lengthWeight = lengthWeight 
        self.rndHubWeight = rndHubWeight  #random Position for radial thread start.
        self.offsetExtSpiral = offsetExtSpiral #offset the radial tread to start from
        self.createCS = createCS  #Create or not the spiral at center
        
        if self.createCS:
            self.centerSpiralThreads = SpiralThread (hub = self.hub, subdiv = self.subdiv, seed= self.seed)
            self.csRotCnt = csRotCnt  #center spiral rotation count
            self.csamax = csamax
            self.csbmax = csbmax
            self.csOffSet = csOffSet
            self.cslgthWei = cslgthWei

        else:
            self.centerSpiralThreads = None
   
        
    ##
    ## @brief      create an array of object used to cast ray on to find anchors points
    ##             Create also an array of GeRayCollider initialise with an object.
    ##
    ## @param      self        
    ## @param      objectList  array of object
    ##
    ## @return     Tru or False if the array could be created
    ##
    def _setObjectList(self, objectList):
        if objectList is None:
            return False
        if len(objectList) == 0:
            return False
        self.objList = objectList
        #creation des RayCollider
        for obj in self.objList:
            if obj is None:
                return False
            self.rcObjList.append(GeRayCollider())
            self.rcObjList[-1].Init(obj)
        return True
        
    ##
    ## @brief      Create an array of radial direction with some randomness (1°-5°)
    ##
    ## @param      self     
    ## @param      raycnt  number of ray to create
    ##
    ## @return     array of vectors
    ##
    def _getRayDirection(self, raycnt):
        #tableau de directions
        #random.seed(self.seed)
        directions = []
        rotation = c4d.utils.Rad(360.0 / float(raycnt))
        for cpt in xrange(raycnt):
            rndangle = random.uniform(0.0,0.5)
            sn, cs = c4d.utils.SinCos(rotation * cpt + c4d.utils.Rad(rndangle))
            #matrice de rotation sur l'axe Z
            rotMatrix = c4d.Matrix( c4d.Vector(0) ,c4d.Vector(cs,sn,0), c4d.Vector(-sn,cs,0) , c4d.Vector(0,0,1) )
            direction  = c4d.Vector(0,1,0) * rotMatrix
            direction.Normalize()
            directions.append( direction)
        return directions

    ##
    ## @brief      cast rays and subcast ray if the ray didn't find any object
    ##
    ## @param      self               
    ## @param      ray_origin_matrix  matrix of the object where the ray should be emmited (will be convert to object space)
    ## @param      distance           length of the ray
    ## @param      raycnt             number of ray to cast
    ## @param      subcast            true or false let know if we are subcasting or not
    ##
    ## @return     array of points if anchors has been found.
    ##
    def _castRay(self, ray_origin_matrix, distance, raycnt, subcast = False):
        anchorPoints = []
        for ray_dir in self._getRayDirection(raycnt):

            for obj, rc in zip(self.objList, self.rcObjList):
                anchorFound = False
                objMg = obj.GetMg()
                ray_origin = ray_origin_matrix.off * ~objMg
                #convertion en coordonnées objet
                obj_ray_dir = (ray_dir * ray_origin_matrix - ray_origin_matrix.off + objMg.off) * ~objMg
                
                rc.Intersect(ray_origin , obj_ray_dir, distance)
                
                if rc.GetIntersectionCount() > 0:
                    if not rc.GetNearestIntersection()["backface"]:
                        anchorPoints.append(rc.GetNearestIntersection()['hitpos'] * objMg * ~self.hub.GetMg())
                        anchorFound = True

            if (not subcast) and not anchorFound:
                    new_rom = c4d.Matrix(ray_origin_matrix)
                    new_rom.off = ray_dir * distance * ray_origin_matrix 
                    
                    subCastTab = None
                    subCastTab = [new_rom.off * ~self.hub.GetMg() ]
                    result = self._castRay(new_rom, distance * 0.1,raycnt,True)
                    
                    subCastTab.append(result)
                    
                    if len(subCastTab[1]) > 0:
                        anchorPoints.append(subCastTab)
                        #print obj


        return anchorPoints
   
    ##
    ## @brief      allow to get the n farrest anchors points without changing the order of the array
    ##
    ## @param      self  
    ## @param      tab   array of anchors points
    ## @param      n     how many number should we keep
    ##
    ## @return     
    ##
    def _sortedValue(self, tab , n = 3):
        #create array of index, value for every value 
        index  = [(i,value) for i, value in enumerate(tab)]
        #sorted the array with the length
        tab_sorted = sorted(index, key = self._getLengthVector)
        #keep only the n bigger
        tab_sorted = tab_sorted[-n:]
        #sort back with the index value (the values get back to their original order)
        tab_sorted = sorted (tab_sorted)
        #rebuild the array with only the value.
        return [v[1] for v in tab_sorted]
        
        
    def _getLongestDistance(self, originPoint = c4d.Vector(0)):
        distToFrame = 0
        ppos = c4d.Vector(0)
        dtoP = 0

        for j in xrange(self.frameThread._getPointCount()):
            ppos = self.frameThread.points[j]
            dtoP = (ppos - originPoint).GetLength()
            if dtoP > distToFrame or distToFrame==0:
                distToFrame = dtoP
        moreLess = float(random.randint(0,3)) / 100.0
        return distToFrame * (0.1 + moreLess)
   
    

    ##
    ## @brief      Compute the spider web so we can display it
    ##
    ## @param      self  
    ##
    ## @return     
    ##
    def _resolveSpiderWeb(self):
        #ancPoints and frame threads
        #
        anchorPoints = []
        hubMg = self.hubMg
        anchorPoints  = self._castRay(self.hubMg,  self.radius, 12)  #default for anchors, 12 rays          
        if len(anchorPoints)  > self.desiredAnchorsCnt:
            #sorted the array and negligate the closer points
            self.frameThread._setPoints( self._sortedValue(anchorPoints, self.desiredAnchorsCnt) ) 
        else:
            self.frameThread._setPoints(anchorPoints)
        if self.frameThread._getAnchorCount() == 0:
            return False

        #support threads
        gv = self._getVector
        anc = self.frameThread._getAnchors()
        ancnt = len(anc)
        ft = [SupportThread(hub = self.hub, subdiv = self.subdiv, seed = self.seed) for i in range (len(anc))]
        
        for i, thr in enumerate(ft):
            if thr._isDone():
                continue
            d = gv(anc[i]).GetLength()
            if d < self.radius * 0.1:
                thr._setDone()
            elif ( gv(anc[i]) - gv(anc[(i+1) % ancnt])).GetLength() < self.radius * 0.35:
                #too short, jump to next anc
                if (gv(anc[(i+1) % ancnt]) - gv(anc[(i+2) % ancnt])).GetLength() < 0:
                    #too short 
                    continue
                else:
                    #Set the link and mark as done, next one also is set as done
                    #end point is 80% with random of 5%
                    ft[i].startPoint = self._getPointPositionAt(
                                            gv(anc[(i-1) % ancnt]),
                                            gv( anc[i]),
                                            self.hubPos,
                                            random.uniform(0.75,0.85))
                    #point d'arrivé dans les 20% avec un petit random de 5%
                    ft[(i) % len(ft)].endPoint = self._getPointPositionAt(  gv(anc[(i+1) % ancnt]),
                                                                            gv(anc[(i+2) % ancnt]),
                                                                            self.hubPos, 
                                                                            random.uniform(0.15,0.25))
                    thr._setDone()
                    ft[(i+1) % len(ft)]._setDone()

            else:
                
                    ft[i].startPoint = self._getPointPositionAt(    gv(anc[(i-1) % ancnt]),
                                                                    gv( anc[i]),
                                                                    self.hubPos,
                                                                    random.uniform(0.75,0.85))
                    ft[(i) % len(ft)].endPoint = self._getPointPositionAt(  gv(anc[(i) % ancnt]),
                                                                            gv(anc[(i+1) % ancnt]),
                                                                            self.hubPos,
                                                                            random.uniform(0.15,0.25))                                                
                
                    thr._setDone()
        self.supportThreads = [t for t in ft if t.startPoint is not None]
        for t in self.supportThreads:
            t._computePoints(c4d.SPLINETYPE_BEZIER)

        ## create the radial thread
        ## 
        ## 

        for i,ray_dir in enumerate(self._getRayDirection(self.radii)):
            ray_end = c4d.Vector(0)
            fnp = c4d.Vector(0) #frame nearest Point
            sfnp = c4d.Vector(0) #support frame nearest Point
            
            distToFrame = 0
            cntsecurity = 0
            ppos = c4d.Vector(0)
            dtoP = 0
            

            while ((distToFrame == 0) or (distToFrame > self._getLongestDistance()* 0.1 )) and (cntsecurity < 10):
                
                ray_end += ray_dir * distToFrame
                
                for j in xrange(self.frameThread._getPointCount()):
                    
                    #ppos = self._getNearestPoint(self.frameThread, ray_end, 1.0 / 10.0)
                    ppos = self.frameThread.points[j]
                    dtoP = (ppos - ray_end).GetLength()
                    if dtoP < distToFrame or distToFrame==0:
                        distToFrame = dtoP
                        fnp = ppos#ray_end
                            
                for st in self.supportThreads:
                    for point in (st.points):
                        ppos = point
                        dtoP = (ppos - ray_end).GetLength()
                        if dtoP < distToFrame or distToFrame==0:
                            distToFrame = dtoP
                            fnp = ppos#ray_end
                
                cntsecurity +=1 

            self.radialThreads.append ( RadialThread(   hub = self.hub,
                                                        subdiv = self.subdiv,
                                                        seed = self.seed,
                                                        startPoint = c4d.Vector(0), 
                                                        endPoint = fnp, 
                                                        ray_dir = (fnp-c4d.Vector(0)).GetNormalized() ,
                                                        length = (fnp-c4d.Vector(0)).GetLength()  
                                                    )    
                                    )    
        for i,ray in enumerate(self.radialThreads):
            if i % 8 == 0:
                #start from center
                ray.startPoint = c4d.Vector(0)
            else:
                #start at random from start of random radial thread created
                ray.startPoint = random.choice(self.radialThreads)._getPoint(random.uniform(0.0,self.rndHubWeight))
            

        ## create the spiral thread
        ## 
        ## 
        distToFrame = 0
        ray_end = 0
        ppos = c4d.Vector(0)
        dtoP = 0
        
        if self.startDistance is not None:
            dFromCenter = self.startDistance
        else:
            dFromCenter = self._getLongestDistance()
        
        self.spiralThreads.points = []
        self.spiralThreads.segments = []

        #equation r = a + b * theta
        a = dFromCenter
        b = 1
        precDir = c4d.Vector(0,1,0)
        rays =  self._getRayDirection(self.radii)
        first = True
        ptcnt = 0

        lenRadialThreads = len(self.radialThreads)
        for tour in xrange(self.rotCnt):
            #for ray in self.radialThreads:
            for i in xrange(lenRadialThreads):
                ray = self.radialThreads[(i+self.offsetExtSpiral) % lenRadialThreads]
                if first:
                    precDir = ray.ray_dir
                    first= False    
                theta = ray.ray_dir.Dot(precDir)
                b = random.uniform(self.bmin,self.bmax) / 100.0  * (1-( (1-a / ray.length)*self.lengthWeight  ) )
                a += random.uniform(self.amin,self.amax) / 100.0 
                r = a + b * theta
                if r <= ray.length:
                    self.spiralThreads.points.append(ray.ray_dir * r)
                    ptcnt +=1
                else:
                    if ptcnt >0:
                        self.spiralThreads.segments.append(ptcnt)
                    ptcnt =0
                precDir = ray.ray_dir


        # create the center spiral 
        #  
        if self.createCS:
            distToFrame = 0
            ray_end = 0
            ppos = c4d.Vector(0)
            dtoP = 0
            dFromCenter = 0.01
            self.centerSpiralThreads.points = []
            self.centerSpiralThreads.segments = []

            #equation r = a + b * theta
            a = dFromCenter
            b = 1
            precDir = c4d.Vector(0,1,0)
            rays =  self._getRayDirection(self.radii)
            first = True
            ptcnt = 0
            
            lenRadialThreads = len(self.radialThreads)
            for tour in xrange(self.csRotCnt):
                #for ray in self.radialThreads:
                for i in xrange(lenRadialThreads):
                    ray = self.radialThreads[(i+self.csOffSet) % lenRadialThreads]
                    if first:
                        precDir = ray.ray_dir
                        first= False    
                    theta = ray.ray_dir.Dot(precDir)
                    
                    b = random.uniform(self.bmin,self.csbmax) / 100.0  * (1-( (1-a / ray.length)*self.cslgthWei  ) )
                    a += random.uniform(self.amin,self.csamax) / 100.0 
                    r = a + b * theta
                    if r <= ray.length:
                        self.centerSpiralThreads.points.append(ray.ray_dir * r)
                        ptcnt +=1
                    else:
                        if ptcnt >0:
                            self.centerSpiralThreads.segments.append(ptcnt)
                        ptcnt =0
                    precDir = ray.ray_dir

class DisplaySpiderWeb(ComputeAndDisplay,object):
    def __init__(self,swg = None, subdiv = 8,**kwargs):
        if swg is None:
            self.swg= None
            return False
        else:
            self.swg = swg
        self.subdiv = subdiv
        super(DisplaySpiderWeb,self).__init__(**kwargs)

    def _getSplines(self):
        container = c4d.BaseObject(c4d.Onull)
        swg = self.swg

        #main frame
        framePoints = swg.frameThread._getPoints()
        frameSpline  = c4d.SplineObject(len(framePoints),c4d.SPLINETYPE_LINEAR)
        #create points
        for i,point in enumerate(framePoints):
            frameSpline.SetPoint(i,point)

        frameSpline[c4d.SPLINEOBJECT_CLOSED] = True
        frameSpline.SetName('Frame Thread')
        frameSpline.InsertUnder(container)
        frameSpline.Message(c4d.MSG_UPDATE)


        #composed Points
        composedPoints = []
        for i in swg.frameThread._getComposedAnchors():
            composedPoints.extend([ [i[0], t] for t in i[1]])
        composedSpline = c4d.SplineObject(len(composedPoints)*2,c4d.SPLINETYPE_LINEAR)
        composedSpline.ResizeObject(len(composedPoints)*2,len(composedPoints))
        for i,comPoint in enumerate(composedPoints):
            composedSpline.SetSegment(i,2,False)
            composedSpline.SetPoint(i+i,comPoint[0])
            composedSpline.SetPoint(i+i+1,comPoint[1])
        composedSpline.SetName('Composed anchor Points')
        composedSpline.InsertUnder(container)
        composedSpline.Message(c4d.MSG_UPDATE)


        #create spline for support threads
        #
        pcnt  = len(swg.supportThreads) * (swg.subdiv + 1 )
        supportThreadsSpline = c4d.SplineObject(pcnt,c4d.SPLINETYPE_LINEAR )
        supportThreadsSpline.ResizeObject(pcnt, len(swg.supportThreads))
        allPoints = []
        for i,spt in enumerate(swg.supportThreads):
            supportThreadsSpline.SetSegment(i, swg.subdiv + 1, False) # x subdivision donne x-1 points
            allPoints.extend(spt._getPoints())
        
        for i, point in enumerate(allPoints):
            supportThreadsSpline.SetPoint(i,point)
        supportThreadsSpline.SetName('support Thread')
        supportThreadsSpline.InsertUnder(container)
        supportThreadsSpline.Message(c4d.MSG_UPDATE)  

        #create radial spline
        #
        radialSpline = c4d.SplineObject(0,c4d.SPLINETYPE_LINEAR)
        radialSpline.ResizeObject(len(swg.radialThreads) * 2,len(swg.radialThreads))

        for i,radthread in enumerate(swg.radialThreads):
            radialSpline.SetSegment(i,2, False)
            radialSpline.SetPoint(i+i,radthread.startPoint)
            radialSpline.SetPoint(i+i+1,radthread.endPoint)
        radialSpline.SetName('Radial Threads')
        radialSpline.InsertUnder(container)
        radialSpline.Message(c4d.MSG_UPDATE)

        ## spiral thread
        ## 
        ## 
        
        spiralThreads = c4d.SplineObject(len(swg.spiralThreads.points),c4d.SPLINETYPE_LINEAR)
        spiralThreads.ResizeObject(len(swg.spiralThreads.points),len(swg.spiralThreads.segments))
        
        for i,poinInSegment in enumerate(swg.spiralThreads.segments):
            spiralThreads.SetSegment(i,poinInSegment,False)
        
        
        for i,p in enumerate(swg.spiralThreads.points):
            spiralThreads.SetPoint(i,p)
        spiralThreads.SetName('Spiral ')
        spiralThreads.InsertUnder(container)
        spiralThreads.Message(c4d.MSG_UPDATE)
        
      
        #return the container      
        container.SetName('Splines')     
        
        ## display center spiral Threads
        if swg.createCS:
            centerSpiralThreads = c4d.SplineObject(len(swg.centerSpiralThreads.points),c4d.SPLINETYPE_LINEAR)
            centerSpiralThreads.ResizeObject(len(swg.centerSpiralThreads.points),len(swg.centerSpiralThreads.segments))
            
            for i,poinInSegment in enumerate(swg.centerSpiralThreads.segments):
                centerSpiralThreads.SetSegment(i,poinInSegment,False)
            
            
            for i,p in enumerate(swg.centerSpiralThreads.points):
                centerSpiralThreads.SetPoint(i,p)
            centerSpiralThreads.SetName('Spiral ')
            centerSpiralThreads.InsertUnder(container)
            centerSpiralThreads.Message(c4d.MSG_UPDATE)

        return container
            


def main():
    start = c4d.GeGetMilliSeconds()

    #container
    container = c4d.BaseObject(c4d.Onull)
    if op[c4d.ID_USERDATA,9]:
        container[c4d.NULLOBJECT_DISPLAY] = 13
    container[c4d.NULLOBJECT_RADIUS] = op[c4d.ID_USERDATA,4] 
    container[c4d.NULLOBJECT_ORIENTATION] = 1
    
    #show plane
    if op[c4d.ID_USERDATA,5]:
        plane = c4d.BaseObject(c4d.Oplane)
        plane[c4d.PRIM_AXIS]=4
        plane[c4d.PRIM_PLANE_SUBW]=1
        plane[c4d.PRIM_PLANE_SUBH]=1
        plane[c4d.PRIM_PLANE_WIDTH] = op[c4d.ID_USERDATA,7]
        plane[c4d.PRIM_PLANE_HEIGHT] = op[c4d.ID_USERDATA,8]
        plane.InsertUnder(container)

             
    swc = ComputeSpiderWeb( hub    = op, 
                            subdiv = op[c4d.ID_USERDATA,11],
                            radius = op[c4d.ID_USERDATA,4],
                            seed   = op[c4d.ID_USERDATA,20],
                            desiredAnchorsCnt = op[c4d.ID_USERDATA,2],
                            radii = op[c4d.ID_USERDATA,3],
                            rotCnt = op[c4d.ID_USERDATA,27],
                            amax=op[c4d.ID_USERDATA,29],
                            bmax=op[c4d.ID_USERDATA,31],
                            startDistance=op[c4d.ID_USERDATA,26],
                            offsetExtSpiral = op[c4d.ID_USERDATA,40],
                            lengthWeight = op[c4d.ID_USERDATA,21],
                            rndHubWeight = op[c4d.ID_USERDATA,13],
                            createCS = op[c4d.ID_USERDATA,39],
                            csRotCnt = op[c4d.ID_USERDATA,25],
                            csOffSet = op[c4d.ID_USERDATA,19],
                            csamax = op[c4d.ID_USERDATA,34],
                            csbmax = op[c4d.ID_USERDATA,37],
                            cslgthWei  = op[c4d.ID_USERDATA,38]
                            
                            

        )
    #get the object list
    objList = []
    for cpt in xrange(op[c4d.ID_USERDATA,1].GetObjectCount()):
        objList.append(op[c4d.ID_USERDATA,1].ObjectFromIndex(doc,cpt))

    if not swc._setObjectList(objList):
        return container
               
  
    ##compute the spider web
    if swc._resolveSpiderWeb() == False:
        return container


    #use the display class to display the result
    displaySpider = DisplaySpiderWeb(swc)
    displaySpider._getSplines().InsertUnder(container) 



    #copy the last hairtag we found on every splines.
    hairtag = op.GetTags()
    hairClone = None
    if hairtag:
        for tag in hairtag:
            if tag.GetType() == c4d.Thairmaterial:
                hairClone =  tag.GetClone()
            
    if hairClone:
        child = container.GetDown().GetDown()
        while child:
                if child.GetType() == c4d.Ospline:
                    child.InsertTag(hairClone.GetClone())
                child = child.GetNext()    

    
    
    end = c4d.GeGetMilliSeconds()
    #print end - start , ' ms to compute'
    #c4d.EventAdd()
    return container

if __name__=='__main__':
    main()
        

