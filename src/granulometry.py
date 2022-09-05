class Granulometry:

    def __init__(self, dm=None, d30=None, d50=None, d90=None, d84tb=None, d84bs=None, Gr=None):
        """
        Reminder : dX is the diameter size such that X percent of the grains have a lower diameter than dX.
        """
        self.dm=dm # Mean diameter
        self.d30=d30
        self.d50=d50 # Median diameter
        self.d90=d90
        self.d84tb=d84tb # TB = travelling bedload
        self.d84bs=d84bs # BS = bed surface
        self.Gr=Gr