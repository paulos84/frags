from django.db import models


class Odor(models.Model):

    CHOICES = (
        ('AL', 'Aldehydic'),
        ('o', 'On loan'),
        ('a', 'Available'),
        ('r', 'Reserved'),
    )



"""
Aldehydic - odor note of the long-chain fatty aldehydes, e.g., fatty-sweaty, ironed laundry, seawater
Animalic  - typical notes from the animal kingdom, e.g., musk, castoreum, skatole, civet, ambergris
Balsamic  - heavy, sweet odors, e.g., cocoa, vanilla, cinnamon, Peru balsam
Camphoraceous  - reminiscent of camphor
Citrus - fresh, stimulating odor of citrus fruits such as lemon or orange
Earthy  - humus-like, reminiscent of humid earth
Fatty  - reminiscent of animal fat and tallow
Floral, flowery generic terms for odors of various flowers
Fruity  - generic term for odors of various fruits
Green  - typical odor of freshly cut grass and leaves
Herbaceous  - noncharacteristic, complex odor of green herbs with, e.g., sage, minty, eucalyptus-like, or earthy nuances
Medicinal  - odor reminiscent of disinfectants, e.g., phenol, lysol, methyl salicylate
Metallic  - typical odor observed near metal surfaces, e.g., brass or steel
Minty  - peppermint-like odor
Mossy  - typical note reminiscent of forests and seaweed
Powdery  - note associated with toilet powders (talcum), diffusively sweet
Resinous  - aromatic odor of tree exudates
Spicy  - generic term for odors of various spices
Waxy  - odor resembling that of candle wax
Woody  - generic term for the odor of wood, e.g., cedarwood, sandalwood
"""
