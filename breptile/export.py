"""STEP export."""

from __future__ import annotations

from OCP.IFSelect import IFSelect_ReturnStatus
from OCP.Message import Message, Message_Alarm
from OCP.Interface import Interface_Static
from OCP.STEPControl import STEPControl_StepModelType, STEPControl_Writer
from OCP.TopoDS import TopoDS_Shape


def write_step(shape: TopoDS_Shape, path: str, schema: str = "AP214") -> None:
    for printer in Message.DefaultMessenger_s().Printers():
        printer.SetTraceLevel(Message_Alarm)  # silence OCC transfer chatter on stdout
    Interface_Static.SetCVal_s("write.step.schema", "AP242DIS" if schema == "AP242" else "AP214DIS")
    writer = STEPControl_Writer()
    status = writer.Transfer(shape, STEPControl_StepModelType.STEPControl_AsIs)
    if status != IFSelect_ReturnStatus.IFSelect_RetDone:
        raise RuntimeError("STEP transfer failed")
    if writer.Write(path) != IFSelect_ReturnStatus.IFSelect_RetDone:
        raise RuntimeError(f"failed to write {path}")
