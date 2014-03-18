import vsp
def geometry_to_clipboard(geometry):
    vsp.CopyGeomToClipboard(geometry.id)

def clipboard_to_geometry(geometry=None):
    if not geometry:
        vsp.PasteGeomClipboard()
    new_children = set(geometry.children) - children
    
    return list(new_children)

class VSP(object):
    @classmethod
    def to_clipboard(cls, geometry):
        #Cache the geometry
        cls._geometry_name = geometry.name
        vsp.CopyGeomToClipboard(geometry)

    @classmethod
    def from_clipboard(cls, parent=None):
        #Pull geometry from cache
        parent_id = ''

        if parent:
            parent_id = parent.id

        vsp.PasteGeomClipboard(parent_id)

        if parent:


        return geometry
        
class Geometry(object):

    def __init__(self, geometry_type, name=None, parent=None):
        if not parent:
            parent = ''

        self._id = vsp.AddGeom(geometry_type, parent.id)

        if name:
            self.name = name
            
        parameter_ids = vsp.GetGeomParmIDs(self._id)

        parameter_types = {
                            0 : float,
                            1 : int,
                            2 : bool,
                            3 : str,
                          }

        for parameter_id in parameter_ids:
            parameter_name = vsp.GetParmName(parameter_id)
            parameter_type = vsp.GetParmType(parameter_id)
            parameter_type = parameter_types[parameter_type]

            setattr(self.__class__, parameter_name.lower(), Parameter(parameter_id, parameter_type))

        self._children = []
        
    def __str__(self):
        string = "{geometry_type}<name={geometry_name}, id={geometry_id}>"
        string = string.format(
                            geometry_type=self.__class__.__name__,
                            geometry_name=self.name,
                            geometry_id=self.id,
                            )
        return string

    @staticmethod
    def to_clipboard(cls):
        
    @property
    def name(self):
        return vsp.GetGeomName(self._id)

    @name.setter
    def name(self, geometry_name):
        vsp.SetGeomName(self._geometry_id, geometry_name)

    @property
    def id(self):
        return self._id

    @property
    def children(self):
        return self._children

    def create_child(self, geometry_type):
        child = Geometry(geometry_type, parent=self)
        self._children.append(child)

        return child

GeometryType = type("GeometryType", (,), {name:value for geometry_type in vsp.GetGeomTypes()})

class Parameter(object):

	def __init__(self, parameter_id, parameter_type):
            self._id = parameter_id
            self._type = parameter_type 

	def __get__(self, instance, owner):
            return vsp.GetParmVal(self._id)

	def __set__(self, instance, value):
            if not isinstance(value, self._type):
                raise RuntimeError("Parameter must be of type {}".format(self._type))
                    
            vsp.SetParmVal(self._id, value, True)
