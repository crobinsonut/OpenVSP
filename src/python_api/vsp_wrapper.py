import vsp

class Geometry(object):

    def __init__(self, geometry_type, name=None, parent=None):
        if not parent:
            parent = ''

        self._id = vsp.AddGeom(geometry_type, parent)

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
        
    def __str__(self):
        string = "{geometry_type}<name={geometry_name}, id={geometry_id}>"
        string = string.format(
                            geometry_type=self.__class__.__name__,
                            geometry_name=self.name,
                            geometry_id=self.id,
                            )
        return string

    @property
    def name(self):
        return vsp.GetGeomName(self._id)

    @name.setter
    def name(self, geometry_name):
        vsp.SetGeomName(self._geometry_id, geometry_name)

    @property
    def id(self):
        return self._id

for geometry_type in vsp.GetGeomTypes():
    setattr(Geometry, geometry_type, geometry_type)

class Pod(Geometry):

    def __init__(self, parent=None, name=None):
        super(Pod, self).__init__(Geometry.POD, name=name, parent=parent)

class Blank(Geometry):

    def __init__(self, name=None, parent=None):
        super(Blank, self).__init__(Geometry.BLANK, name=name, parent=parent)

class Fuselage(Geometry):

    def __init__(self, name=None, parent=None):
        super(Fuselage, self).__init__(Geometry.FUSELAGE, name=name, parent=parent)

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
