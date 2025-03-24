from abc import ABC, abstractmethod


class SpectralOperatorStrategy(ABC):

    @abstractmethod
    def GetMatrix(self, energy):
        pass
        

    @abstractmethod
    def SpectralFunction(self, x, energy):
        pass
        
